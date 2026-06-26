#========================================================================================#
"""
	NeutralDrift

This model runs a very simple simulation called a Moran process. At each step, individuals change
their state to that of a random neighbour. 

Author: Niall Palfreyman, March 2025.
"""
module NeutralDrift

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Moran

A Moran agent has an attribute (in this case `altruist`) that is determined at each step by the
distribution of attributes of its neighbours.
"""
@agent struct Moran(ContinuousAgent{2,Float64})
	altruist::Bool							# Type of the agent
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	neutraldrift( kwargs)

Create and initialise the NeutralDrift model.
"""
function neutraldrift(;
	altruism_cost = 0.0,					# Individual cost of acting altruistically
)
	extent = (30,30)
	properties = Dict(
		:altruism_cost => altruism_cost,
	)

	moranworld = StandardABM( Moran, ContinuousSpace(extent; spacing = 1.0);
		agent_step!,
		properties
	)

	for idx in CartesianIndices(extent)
		add_agent!( idx, moranworld, (0,0), rand(Bool))
	end

	return moranworld
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( moran, world)

On each step, each Moran agent replaces its altruistic attribute by that of a random neighbour.
"""
function agent_step!( moran::Moran, world::ABM)
	nbr = random_nearby_agent( moran, world, 1.5)
	moran.altruist =
		(nbr.altruist && rand() < 1 - world.altruism_cost) ? true : false
end

#-----------------------------------------------------------------------------------------
"""
	acolour( moran::Moran)

Select the agent's colour on the basis of their altruism.
"""
agcolour( moran::Moran) = moran.altruist ? :lime : :blue

#-----------------------------------------------------------------------------------------
"""
	demo()

Run a simulation of the NeutralDrift model.
"""
function demo()
	params = Dict(
		:altruism_cost => 0:1e-5:1e-3
	)
	plotkwargs = (
        # To-do: Specify plotting keyword arguments
		agent_color = agcolour,
		agent_size = 20,
	)

	playground, = abmplayground( neutraldrift; params, plotkwargs...)

	playground
end

end