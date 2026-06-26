#========================================================================================#
"""
	Stabilisation

Demonstrate the stabilisation of movement within a dynamical system.

Author: Niall Palfreyman, March 2025.
"""
module Stabilisation

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Particle

A Particle simply moves around with a speed, however this speed is influenced by the presence of
other nearby Particles.
"""
@agent struct Particle(ContinuousAgent{2,Float64})
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	stabilisation()

Create the stabilisation model.
"""
function stabilisation(;
	binding_radius = 0.2,			# Radius of attraction between Particles
)
	extent = (20,20)
	properties = Dict(
		:binding_radius => binding_radius,
	)

	world = StandardABM( Particle, ContinuousSpace(extent, spacing=0.5);
		agent_step!, properties
	)
	
	for _ in 1:prod(extent)
		# Random facing direction:
		theta = 2pi*rand()
		add_agent!( world, (cos(theta),sin(theta)))
	end

	return world
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( particle, world)

The particle moves in its facing direction with a speed that adjusts according to the number of
nearby particles, and then wiggles to face in a new random direction.
"""
function agent_step!( particle, world)
	# Chance of moving is lower if other particles are nearby:
	if rand() < (1 / (1 + length(collect(nearby_agents(particle,world,world.binding_radius)))^2))
		move_agent!(particle, world)
	end
	wiggle!(particle,pi/9)
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Set up a playground for a simple stabilisation process.
"""
function demo()
	params = Dict(
		:binding_radius => 0.05:0.05:1,
	)

	plotkwargs = (
		agent_color=:green,
		agent_marker=(ag->wedge(ag,0.5)),
	)

	playground, = abmplayground( stabilisation; params, plotkwargs...)

	playground
end

end