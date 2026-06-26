"""
	Suboptimisation

Demonstrate the notion of suboptimisation using two objective functions: a set of three valleys
and a member of the De Jong test suite.

Author: Niall Palfreyman, March 2025.
"""
module Suboptimisation

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Scout

Scouts fly around the world seeking to minimise the value of an objective function.
"""
@agent struct Scout(ContinuousAgent{2,Float64})
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
Initialise the Suboptimisation model, also setting up an objective function which is either
easy (valleys) or difficult (De Jong 2).
"""
function suboptimisation(;
	difficult::Bool=false,
	speed::Float64=0.1,
)
	extent = (80,80)	# Extent of world
	pPop = 0.1			# Proportion of populated locations

	properties = Dict(
		:objective => (difficult ? dejong2(extent) : valleys(extent)),
		:pPop => pPop,
		:difficult => difficult,
		:speed => speed,
	)
	
	world = StandardABM( Scout, ContinuousSpace(extent, spacing=1.0);
		agent_step!,
		properties
	)
	
	for _ in 1:pPop*prod(extent)
		theta  = 2pi*rand()
		add_agent!( world, (cos(theta),sin(theta)))
	end

	world
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!(scout,world)

Scout takes a small step in the direction of a nearby location with a smaller objective value than
at the Scout's own location.
"""
function agent_step!( scout, world)
	scout.vel = .-gradient(scout.pos,world.objective,world)
	move_agent!(scout,world,world.speed)
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Demonstrate the Suboptimisation model.
"""
function demo()
	params = Dict(
		:difficult => false:true,
		:speed => 0:0.1:5,
	)
	plotkwargs = (
		agent_marker = wedge,
		agent_color = :red,
		add_colorbar=false,
		heatarray=:objective,
	)

	playground, = abmplayground( suboptimisation; params, plotkwargs...)

	playground
end

end