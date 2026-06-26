#========================================================================================#
"""
	MultiLevel

This model demonstrates how agents can recruit each other into a communal group that cooperatively
solves a communal problem. In evolution, we call this recruitment in the interest of survival
Multi-Level Selection; in studying behaviour, we call it Cooperation.

This model is adapted from Wilensky's EACH model of cooperation at:
	http://ccl.northwestern.edu/netlogo/models/Altruism

Author: Niall Palfreyman, March 2025.
"""
module MultiLevel

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Individual

An Individual agent has the trait `cooperative`, whose Bool value is genetically inherited, and
the trait `fitness`, which is calculated from its interaction with neighbours.
"""
@agent struct Individual(ContinuousAgent{2,Float64})
	cooperative::Bool									# Is the Individual cooperative?
	fitness::Float64									# Individual's fitness
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	multi_level( kwargs)

Create and initialise a MultiLevel selection model.
"""
function multi_level(;
	cost_cooperativity = 0.1,					# Cost to the individual of cooperating
	bene_cooperativity = 0.6,					# Benefit of having a cooperating neighbour
	prob_dying = 0.1,							# What is the probability of an individual dying?
	prob_birth = 1.0,							# Probability of being able to give birth
)
	extent = (51,51)
	prob_populated = 0.5
	prob_cooperativity = 0.5

	properties = Dict(
		:cost_cooperativity => cost_cooperativity,
		:bene_cooperativity => bene_cooperativity,
		:prob_dying => prob_dying,
		:prob_birth => prob_birth,
		:resource => zeros(Float64,extent),
	)

	ml = StandardABM( Individual, ContinuousSpace(extent; spacing = 1.0);
		properties, model_step!, agent_step!
	)

	individuals = rand(round(Int,prob_populated*prod(extent))) .< prob_cooperativity
	map(individuals) do cooperative
		add_agent!( ml, (0,0), cooperative, 0.0)
	end

	return ml
end

#-----------------------------------------------------------------------------------------
"""
	model_step!( ml)

Use individuals' fitness to decide whether they give birth or die.
"""
function model_step!( ml::ABM)
	for agent in collect(allagents(ml))
		nbrs = collect(nearby_agents(agent,ml,1.5))
		if length(nbrs) < 8 && rand() < ml.prob_birth
			# Give birth to new neighbour with locally prevalent cooperativity trait:
			nbhd_fitness = agent.fitness
			coop_fitness = agent.cooperative ? agent.fitness : 0.0

			for nbr in nbrs
				nbhd_fitness += nbr.fitness				# Total fitness available in neighbourhood
				if nbr.cooperative
					coop_fitness += nbr.fitness			# Total fitness of cooperators in nbhood
				end
			end																	# Baby inherits from
			cooperate = rand() < coop_fitness/nbhd_fitness						# fittest nbhd group
			birthpos = normalize_position(agent.pos.+Tuple(2rand(2).-1),ml)		# Give birth nearby

			add_agent!( birthpos, ml, (0,0), cooperate, 0.0)					# A child is born!
		end

		if rand() < ml.prob_dying || rand() > agent.fitness						# Fitness reduces
			remove_agent!(agent,ml)												# death rate
		end
	end
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( inidvidual, ml)

Calculate individual's fitness from the number of cooperators in its Moore 8-neighbourhood.
"""
function agent_step!( individual::Individual, ml::ABM)
	nbrs = nearby_agents(individual,ml,1.5)									# Moore neighbourhood
	nbhd_sz = 8																# contains 8 neighbours
	n_cooperators = count(map(agent->agent.cooperative,nbrs))				# Gain fitness from each
	individual.fitness = 0.1 + ml.bene_cooperativity*n_cooperators/nbhd_sz	# cooperative neighbour
																			# ... and also from myself
	if individual.cooperative												# if I am cooperative, but
		individual.fitness += ml.bene_cooperativity/nbhd_sz-ml.cost_cooperativity	# I pay the cost!
	end
end

#-----------------------------------------------------------------------------------------
"""
	acolour( individual::Individual)

Select the individual's colour on the basis of its cooperativity.
"""
acolour( individual::Individual) = individual.cooperative ? :lime : :blue

#-----------------------------------------------------------------------------------------
"""
	demo()

Run a simulation of the MultiLevel model.
"""
function demo()
	params = Dict(
		:cost_cooperativity => 0:0.01:1,
		:bene_cooperativity => 0:0.01:1,
		:prob_dying => 0:0.01:0.2,
		:prob_birth => 0.8:0.01:1,
	)
	plotkwargs = (
		agent_color = acolour,
		agent_size = 20,
		adata=[(ag->ag.cooperative,count),(ag->!ag.cooperative,count)],
		alabels=["Cooperators","Defectors"],
	)

	playground, = abmplayground( multi_level; params, plotkwargs...)

	playground
end

end