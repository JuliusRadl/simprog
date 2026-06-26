#========================================================================================#
"""
	Genetic structure

A simulation to investigate the use of genetic structure in agents to help them solve problems.
	
Author: Niall Palfreyman, March 2025.
"""
module GeneticStructure

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, Random, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Turtle

Turtles sit in one place, unmoving, and are selected for reproduction depending on how well
they minimise a particular objective function: dissonance.
"""
@agent struct Turtle(ContinuousAgent{2,Float64})
	genome::Vector{Char}			# Turtle's character genome
	dissonance::Float64				# Dissonance of the Turtle with its environment
	colour::Symbol					# Colour indicating turtle's dissonance status
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	genetic_structure()

Initialise the Genetic Structure simulation.
"""
function genetic_structure(;
	crossover = true,
	mu_rate = 1e-3,
	elitism = -1,
)
	extent = (40,40)								# Extent of the world
#	nPop = round(Int,0.3*prod(extent))				# Population size
	nPop = 10
	properties = Dict(
		:crossover		=> crossover,				# Are we using crossover?
		:mu_rate		=> mu_rate,					# Probability of mutating a locus
		:elitism		=> elitism,					# How much dissonance are we eliminating?
		:nPop			=> nPop,					# (Constant) population size
		:roulette		=> cumsum(ones(nPop)/nPop),	# Roulette wheel selector
		:target			=> Char.([
			84,104,101,32,115,116,97,114,114,121,32,104,101,97,118,101,110,115,32,
			97,98,111,118,101,32,109,101,44,32,97,110,100,32,116,104,101,32,
			109,111,114,97,108,32,108,97,119,32,119,105,116,104,105,110,32,109,101
		]),
		:alphabet		=> Char.(32:122),			# Collection of available alleles
		:mean			=> 1.0,						# Mean dissonance of population
		:minID			=> 1,						# ID of turtle with minimum dissonance
	)

	gs = StandardABM( Turtle, ContinuousSpace(extent);
		model_step!, properties
	)
	genome_length = length(gs.target)
	for _ in 1:nPop
		lilith = add_agent!( gs, (1,1),				# Ur-Turtle has ...
			rand( gs.alphabet, genome_length),		# Genome (candidate Vector{Char})
			1.0, :white								# Provisional dissonance and fitness type
		)
		set_dissonance!( lilith, gs)					# Set lilith's dissonance/colour
	end

	gs
end

#-----------------------------------------------------------------------------------------
"""
	model_step!( gs)

Perform dissonance-based selection on the Turtle population to decide who gets children.
"""
function model_step!( gs::ABM)
	allturtles = collect(allagents(gs))
	dissonances = (ag->ag.dissonance).(allturtles)
	ids = (ag->ag.id).(allturtles)
	gs.mean = mean(dissonances)
	gs.minID = allturtles[findmin(dissonances)[2]].id

	sigma = std(dissonances)					# Standard deviation
	if sigma == 0
		# Singular case: all evaluations were equal to the mean:
		fitnesses = ones(gs.nPop)
	else
		# Otherwise perform sigma-scaling and exclude lowest elitism percentile of population:
		fitnesses = (gs.mean .- dissonances)/sigma .- gs.elitism
		fitnesses[fitnesses .<= 0] .= 0
	end
	roulette = cumsum( fitnesses)/sum(fitnesses)		# Roulette wheel of dissonance-biased slots
	throws = rem.( rand() .+ (1:gs.nPop)./gs.nPop, 1)	# Statistical Uniformly Sampled ball throws

	# Throw ball onto roulette wheel to create randomised list of fitness-selected parent ids:
	parents = zeros(Integer,gs.nPop)
	for parent in 1:gs.nPop
		for slot in 1:gs.nPop
			if throws[parent] <= roulette[slot]
				parents[parent] = ids[slot]
				break
			end
		end
	end
	shuffle!(parents)

	# 1st half of parents are Mummies; 2nd half are Daddies:
	nMatings = gs.nPop ÷ 2
	for mating in 1:nMatings
		reproduce( parents[mating], parents[nMatings+mating], gs)
	end
end

#-----------------------------------------------------------------------------------------
"""
	reproduce( mum, dad, gs)

Give birth to a baby with similar attributes to mum and dad, but possibly genetically mutated with
probability gs.mu_rate. NOTE: This is the core of the genetic algorithm!
"""
function reproduce( mum, dad, gs)
	mummy = gs[mum]
	daddy = gs[dad]
	# Start with crossover ...
	len = length(mummy.genome)
	xpt = gs.crossover ? rand(1:len-1) : len					# Crossover point
	billy = deepcopy(vcat(daddy.genome[1:xpt],mummy.genome[xpt+1:end]))	# Baby Billy's genome
	sally = deepcopy(vcat(mummy.genome[1:xpt],daddy.genome[xpt+1:end]))	# Baby Sally's genome

	# ... then mutate ...
	mutated_loci = rand(len).<gs.mu_rate
	billy[mutated_loci] = rand(gs.alphabet,count(mutated_loci))
	mutated_loci = rand(len).<gs.mu_rate
	sally[mutated_loci] = rand(gs.alphabet,count(mutated_loci))

	# ... transfer children Billy and Sally to replace parents:
	mummy.genome[:] = sally[:]; set_dissonance!( mummy, gs)
	daddy.genome[:] = billy[:]; set_dissonance!( daddy, gs)
end

#-----------------------------------------------------------------------------------------
"""
	set_dissonance!( turtle, gs)

Set turtle's initial colour to reflect its dissonance with respect to the gs model's
target string. A turtle's dissonance is the mean Hamming distance per character of its genome from
the target.
"""
function set_dissonance!( turtle::Turtle, gs)
	turtle.dissonance = dissonance( turtle.genome, gs.target)	# Mean Manhattan difference per locus

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

	playground,abmplt = abmplayground( genetic_structure; params, plotkwargs...)

	best_string = lift( (wld->String(wld[wld.minID].genome)), abmplt.model)
	text!( 40,30, text=best_string, color=:red, fontsize=16, align=(:left,:top))

	playground
end

end