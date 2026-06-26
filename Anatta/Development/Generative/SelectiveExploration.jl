#========================================================================================#
"""
	Selective exploration

A simulation to investigate the use of genetic selection to help agents solve problems. Unlike
SelectiveUnits in the Selective Search simulation, ExploratoryAgents in this simulation possess
Agency. That is, they move actively around the world, Exploring their environment and gaining
energy from their interactions with other agents.
	
Author: Niall Palfreyman, January 2026.
"""
module SelectiveExploration

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	ExploratoryAgent

ExploratoryAgents chase around the world looking for a partner to reproduce with. They only survive
as long as they have positive energy, and their initial energy depends on how well they minimise a
particular objective function: dissonance.
"""
@agent struct ExploratoryAgent(ContinuousAgent{2,Float64})
	genome::Vector{Char}			# ExploratoryAgent's character genome
	energy::Float64					# Current energy of the ExploratoryAgent
	dissonance::Float64				# Dissonance of the ExploratoryAgent with its environment
	colour::Symbol					# Colour indicating ExploratoryAgent's energy status
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	selective_exploration()

Initialise the Selective Exploration simulation.
"""
function selective_exploration(;
	crossover = true,
	mu_rate = 1e-3,
	unimodal_fitness = true,
	difficulty = 0.0,
	benefit = 0.9,										# Benefit of lower dissonance (relative to cost of living
)
	pPop = 0.3											# Probability of a cell becoming populated
	extent = (40,40)									# Extent of the world
	properties = Dict(
		:crossover				=> crossover,			# Are we using crossover?
		:mu_rate				=> mu_rate,				# Probability of mutating a locus
		:unimodal_fitness		=> unimodal_fitness,	# Use unimodal fitness landscape?
		:difficulty				=> difficulty,			# Difficulty of search problem
		:benefit				=> benefit,				# Benefit of lower dissonance (relative to cost of living)
		:life_energy			=> 10.0,				# Average initial energy of ExploratoryAgents
		:living_cost			=> 1.0,					# How much energy does living cost?
		:birth_cost				=> 1.0,					# How much energy does giving birth cost?
		:search_speed			=> 1.0,					# How fast should ExploratoryAgents move?
		:interaction_distance	=> 5.0,					# From what distance can ExploratoryAgents interact?
		:target					=> quotation,			# Target string to be evolved
		:alphabet				=> quot_alphabet,		# Collection of available alleles
		:max_pop 				=> pPop*prod(extent),	# Maximum (approx) allowed population
		:minID					=> 1,					# ID of ExploratoryAgent with minimum dissonance
		:mean					=> 0.0,					# Mean dissonance of population
	)
	se = StandardABM( ExploratoryAgent, ContinuousSpace(extent);
		model_step!, agent_step!, properties
	)
	genome_length = length(se.target)
	for _ in 1:se.max_pop
		θ = 2π * rand()
		earnie = add_agent!( se, (cos(θ),sin(θ)),	# Ur-ExploratoryAgent has ...
			rand( se.alphabet, genome_length),		# Genome (candidate Vector{Char})
			2rand() * se.life_energy,				# Initial energy
			1.0,									# Provisional dissonance
			:white,									# Provisional colour
		)
		set_dissonance!( earnie, se)				# Set earnie's dissonance/colour
	end

	se
end

#-----------------------------------------------------------------------------------------
"""
	model_step!( se)

Locate agent with current minimum dissonance and store its id in se.minID.
"""
function model_step!( se::ABM)
	allExploratoryAgents = collect(allagents(se))
	dissonances = (ag->ag.dissonance).(allExploratoryAgents)
	se.minID = allExploratoryAgents[findmin(dissonances)[2]].id
	se.mean = mean(dissonances)
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( ExploratoryAgent, se)

One step in the life of a ExploratoryAgent with senetic structure.
"""
function agent_step!( earnie::ExploratoryAgent, se)
	walk!( earnie, se)									# Move around the world
	explore!( earnie, se)								# Gain energy from interaction with others
	if earnie.energy < 0								# Die, possibly
		remove_agent!(earnie, se)
		return
	end
	reproduce!( earnie, se)								# Reproduce if possible
end

#-----------------------------------------------------------------------------------------
"""
	walk!( ExploratoryAgent, se)

Walk approximately forwards with speed up to the maximum, and age by one cost of living unit.
"""
function walk!( earnie::ExploratoryAgent, se)
	wiggle!(earnie)
	move_agent!(earnie,se,rand()*se.search_speed)
	earnie.energy -= se.living_cost						# Cost of living
end

#-----------------------------------------------------------------------------------------
"""
	explore!( ExploratoryAgent, se)

Interact with other nearby agents, possibly profiting from this interaction.
"""
function explore!( earnie::ExploratoryAgent, se)
	nbr = random_nearby_agent(earnie,se,se.interaction_distance)
	if nbr !== nothing &&
			spiky(1-earnie.dissonance,se.difficulty) > spiky(1-nbr.dissonance,se.difficulty)
		# Reward earnie for having locally lower dissonance:
		earnie.energy += se.living_cost * se.benefit
	end
end

#-----------------------------------------------------------------------------------------
"""
	reproduce!( ExploratoryAgent, se)

Give birth to a baby with similar attributes to self, but possibly genetically mutated with
probability se.mu_rate. NOTE: This is the core of the genetic algorithm!
"""
function reproduce!( mummy::ExploratoryAgent, se)
	if nagents(se) < se.max_pop && mummy.energy > se.birth_cost
		# I've got enough energy and there's sufficient space to have kids:
		daddy = random_nearby_agent(mummy,se,se.interaction_distance)
		if daddy !== nothing && daddy.energy > se.birth_cost
			# Let's mate! Start with crossover ...
			len = length(mummy.genome)
			xpt = se.crossover ? rand(1:len-1) : len					# Crossover point
			billy = vcat(daddy.genome[1:xpt],mummy.genome[xpt+1:end])	# Baby Billy's genome
			sally = vcat(mummy.genome[1:xpt],daddy.genome[xpt+1:end])	# Baby Sally's genome

			# ... then mutate ...
			mutated_loci = rand(len).<se.mu_rate
			billy[mutated_loci] = rand(se.alphabet,count(mutated_loci))
			mutated_loci = rand(len).<se.mu_rate
			sally[mutated_loci] = rand(se.alphabet,count(mutated_loci))

			# ... construct children Billy and Sally ...
			θ,ϕ = 2π*rand(2)
			set_dissonance!(
				add_agent!( se, (cos(θ),sin(θ)), billy, 2rand() * se.life_energy, 1.0, :white),
				se
			)
			set_dissonance!(
				add_agent!( se, (cos(ϕ),sin(ϕ)), sally, 2rand() * se.life_energy, 1.0, :white),
				se
			)

			# ... and deplete parents:
			mummy.energy -= se.birth_cost
			daddy.energy -= se.birth_cost
		end
	end
end

#-----------------------------------------------------------------------------------------
"""
	set_dissonance!( ExploratoryAgent, se)

Set ExploratoryAgent's dissonance and its colour to reflect its dissonance with respect to the se model's
target string. A ExploratoryAgent's dissonance is the mean Hamming distance per character of its genome from
the target.
"""
function set_dissonance!( earnie::ExploratoryAgent, se)
	if se.unimodal_fitness
		diss = dissonance( earnie.genome, se)			# Mean Manhattan difference per locus
	else
		diss = bimodal_dissonance( earnie.genome, se)	# Bimodal dissonance
	end
	earnie.dissonance = min(1.0, diss)

	# Colour ExploratoryAgent to indicate dissonance level
	if earnie.dissonance > 0.5
		earnie.colour = :black
	elseif earnie.dissonance > 0.2
		earnie.colour = :blue
	elseif earnie.dissonance > 0.02
		earnie.colour = :red
	elseif earnie.dissonance > 0.002
		earnie.colour = :orange
	else
		earnie.colour = :yellow
	end
end

#-----------------------------------------------------------------------------------------
"""
	dissonance( candidate::Vector{Char}, se)

The dissonance of a candidate string with respect to a target string (of the same length as the
candidate) is its Hamming distance from the target, divided by the length of the target - so, the
mean Hamming distance per string character.
"""
function dissonance( candidate::Vector{Char}, se)
	divisor = length(se.target) * (length(se.alphabet) - 1)
	sum(abs.(candidate .- se.target))/divisor
end

#-----------------------------------------------------------------------------------------
"""
	bimodal_dissonance( candidate::Vector{Char}, se)

This dissonance function still has the same global minimum at 0 as dissonance() when the
candidate and target strings are equal to each other, but introduces a second local minimum.
"""
function bimodal_dissonance( candidate::Vector{Char}, se)
	diss = dissonance(candidate,se)
	diss + 0.2sin(4pi*diss)
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Demonstrate a simple genetic algorithm.
"""
function demo()
	params = Dict(
		:mu_rate => 0:0.001:0.1,
		:benefit => 0.0:0.1:1.5,
		:unimodal_fitness => (false, true),
		:difficulty => 0:.01:1,
	)
	plotkwargs = (
		dt=1:200,
		agent_color=(ExploratoryAgent->ExploratoryAgent.colour),
		agent_size=20,
		mdata = [(m->m.mean), (m->m[m.minID].dissonance)],
		mlabels = ["Mean dissonance", "Minimum dissonance"],
	)

	playground,abmplt = abmplayground( selective_exploration; params, plotkwargs...)

	best_string = lift( (wld->String(wld[wld.minID].genome)), abmplt.model)
	text!( 10,1.0, text=best_string, color=:red, fontsize=16, align=(:left,:top))

	playground
end

end