#========================================================================================#
"""
	Genetic exploration

A simulation to investigate the use of genetically constrained exploration in agents to help them
solve problems.
	
Author: Niall Palfreyman, March 2025.
"""
module GeneticExploration

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Turtle

Turtles chase around the world looking for a partner to reproduce with. They only survive as long
as they have positive energy, and their initial energy depends on how well they minimise a
particular objective function: dissonance.
"""
@agent struct Turtle(ContinuousAgent{2,Float64})
	genome::Vector{Char}			# Turtle's character genome
	energy::Float64					# Current energy of the Turtle
	dissonance::Float64				# Dissonance of the Turtle with its environment
	colour::Symbol					# Colour indicating turtle's energy status
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	genetic_exploration()

Initialise the Genetic Structure simulation.
"""
function genetic_exploration(;
	crossover = true,
	mu_rate = 1e-3,
	elitism = -1,
)
	pPop = 0.2										# Probability of a cell becoming populated
	extent = (40,40)								# Extent of the world
	properties = Dict(
		:crossover		=> crossover,				# Are we using crossover?
		:mu_rate		=> mu_rate,					# Probability of mutating a locus
		:elitism		=> elitism,					# How much dissonance are we eliminating?
		:life_energy	=> 10.0,					# Average initial energy of turtles
		:living_cost	=> 0.1,						# How much of life energy does living cost?
		:search_speed	=> 1.0,						# How fast should turtles move?
		:target			=> Char.([
			84,104,101,32,115,116,97,114,114,121,32,104,101,97,118,101,110,115,32,
			97,98,111,118,101,32,109,101,44,32,97,110,100,32,116,104,101,32,
			109,111,114,97,108,32,108,97,119,32,119,105,116,104,105,110,32,109,101
		]),
		:alphabet		=> Char.(32:122),			# Collection of available alleles
		:max_pop 		=> pPop*prod(extent),		# Maximum (approx) allowed population
		:mean			=> 1.0,						# Mean dissonance of population
		:sigma			=> 0.0,						# Standard deviation of dissonance
		:minID			=> 1,						# ID of turtle with minimum dissonance
	)

	ge = StandardABM( Turtle, ContinuousSpace(extent);
		model_step!, agent_step!, properties
	)
	genome_length = length(ge.target)
	for _ in 1:ge.max_pop/2
		θ = 2π * rand()
		lilith = add_agent!( ge, (cos(θ),sin(θ)),	# Ur-Turtle has ...
			rand( ge.alphabet, genome_length),		# Genome (candidate Vector{Char})
			1.0, 1.0,								# Provisional energy and dissonance
			:white,									# Provisional fitness type
		)
		set_fitness!( lilith, ge)					# Set lilith's fitness/colour
	end

	ge
end

#-----------------------------------------------------------------------------------------
"""
	model_step!( ge)

Calculate current mean and standard deviation dissonance to enable fitness calculations.
"""
function model_step!( ge::ABM)
	allturtles = collect(allagents(ge))
	dissonances = (ag->ag.dissonance).(allturtles)
	ge.sigma = std(dissonances)
	ge.mean = mean(dissonances)
	ge.minID = allturtles[findmin(dissonances)[2]].id
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( turtle, ge)

One step in the life of a Turtle with genetic structure.
"""
function agent_step!( turtle::Turtle, ge)
	walk!( turtle, ge)
	if turtle.energy < 0
		remove_agent!(turtle, ge)
		return
	end
	reproduce( turtle, ge)
end

#-----------------------------------------------------------------------------------------
"""
	walk!( turtle, ge)

Walk approximately forwards with speed up to the maximum, and age by one cost of living unit.
"""
function walk!( turtle::Turtle, ge)
	wiggle!(turtle)
	move_agent!(turtle,ge,rand()*ge.search_speed)
	turtle.energy -= ge.living_cost * ge.life_energy
end

#-----------------------------------------------------------------------------------------
"""
	reproduce( turtle, ge)

Give birth to a baby with similar attributes to self, but possibly genetically mutated with
probability ge.mu_rate. NOTE: This is the core of the genetic algorithm!
"""
function reproduce( mummy::Turtle, ge)
	birth_cost = ge.living_cost * ge.life_energy
	if nagents(ge) < ge.max_pop && mummy.energy > birth_cost
		# I've got enough energy and there's sufficient space to have kids:
		daddy = random_nearby_agent(mummy,ge)
		if daddy !== nothing && daddy.energy > birth_cost
			# Let's mate! Start with crossover ...
			len = length(mummy.genome)
			xpt = ge.crossover ? rand(1:len-1) : len					# Crossover point
			billy = vcat(daddy.genome[1:xpt],mummy.genome[xpt+1:end])	# Baby Billy's genome
			sally = vcat(mummy.genome[1:xpt],daddy.genome[xpt+1:end])	# Baby Sally's genome

			# ... then mutate ...
			mutated_loci = rand(len).<ge.mu_rate
			billy[mutated_loci] = rand(ge.alphabet,count(mutated_loci))
			mutated_loci = rand(len).<ge.mu_rate
			sally[mutated_loci] = rand(ge.alphabet,count(mutated_loci))

			# ... construct children Billy and Sally ...
			θ,ϕ = 2π*rand(2)
			set_fitness!( add_agent!( ge, (cos(θ),sin(θ)), billy, 1.0, 1.0, :white), ge)
			set_fitness!( add_agent!( ge, (cos(ϕ),sin(ϕ)), sally, 1.0, 1.0, :white), ge)

			# ... and deplete parents:
			mummy.energy -= birth_cost
			daddy.energy -= birth_cost
		end
	end
end

#-----------------------------------------------------------------------------------------
"""
	set_fitness!( turtle, ge)

Set turtle's initial energy and its colour to reflect its dissonance with respect to the ge model's
target string. A turtle's dissonance is the mean Hamming distance per character of its genome from
the target.

Note that an individual's fitness is context-dependent: it is calculated from the individual's
standing with respect to the entire current population (summarised by ge.mean and ge.sigma).
"""
function set_fitness!( turtle::Turtle, ge)
	turtle.dissonance = dissonance( turtle.genome, ge.target)	# Mean Manhattan difference per locus
	if ge.sigma==0
		# Really bizarre case at beginning of run - all turtles are equally dissonant:
		fitness = 1.0
	else
		# Otherwise perform sigma-scaling and exclude lowest percentile of the population:
		fitness = (ge.mean - turtle.dissonance)/(ge.sigma) - ge.elitism
	end
	turtle.energy = fitness>0 ? (fitness*ge.life_energy) : -1

	# Colour turtle to indicate dissonance level
	if turtle.dissonance > 5.0
		turtle.colour = :black
	elseif turtle.dissonance > 1.0
		turtle.colour = :blue
	elseif turtle.dissonance > 0.5
		turtle.colour = :red
	elseif turtle.dissonance > 0.2
		turtle.colour = :orange
	else
		turtle.colour = :yellow
	end
end

#-----------------------------------------------------------------------------------------
"""
	dissonance( candidate::Vector{Char}, target::Vector{Char})

The dissonance of a candidate string with respect to a target string (of the same length as the
candidate) is its Hamming distance from the target, divided by the length of the target - so, the
mean Hamming distance per string character.
"""
function dissonance( candidate::Vector{Char}, target::Vector{Char})
	sum(abs.(candidate - target))/length(target)
end

#-----------------------------------------------------------------------------------------
"""
	local_dissonance( candidate::Vector{Char}, target::Vector{Char})

This dissonance function still has the same global minimum at 0 as dissonance() when the
candidate and target strings are equal to each other.
"""
function local_dissonance( candidate::Vector{Char}, target::Vector{Char})
	diss = dissonance(candidate,target)
	0.2diss^3 - 3.9diss^2 + 20diss
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Demonstrate a simple genetic algorithm.
"""
function demo()
	params = Dict(
		:crossover => false:true,
		:mu_rate => 0:0.001:0.1,
		:elitism => -3:3,
	)
	plotkwargs = (
		agent_color=(turtle->turtle.colour),
		agent_size=20,
		mdata = [(m->m.mean), (m->m[m.minID].dissonance)],
		mlabels = ["Mean dissonance", "Minimum dissonance"],
	)

	playground,abmplt = abmplayground( genetic_exploration; params, plotkwargs...)

	best_string = lift( (wld->String(wld[wld.minID].genome)), abmplt.model)
	text!( 40,30, text=best_string, color=:red, fontsize=16, align=(:left,:top))

	playground
end

end