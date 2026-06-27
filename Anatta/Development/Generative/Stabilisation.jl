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
	sticky_x = 0,
	sticky_y = 0
)
	extent = (20, 20)
	properties = Dict(
		:binding_radius => binding_radius,
		:sticky_x => sticky_x,
		:sticky_y => sticky_y,
		:extent => extent
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
		# chance of moving is lower if near sticky location
		# distance between agent pos and sticky pos
		x = world.sticky_x
		y = world.sticky_y
		sticky_pos = (x, y)
		dist = Agents.euclidean_distance(particle.pos, sticky_pos, world)

		# distance as proportion of maximum distance possible in the space
		# (distance between the center and a corner, since the space is periodic)
		max_dist = ((world.extent[1]/2)^2 + (world.extent[2]/2)^2)^0.5
		
		# relative distance squared (exponentially stronger effect when closer)
		proximity = (1 - (dist / max_dist)^2)
		if rand() > proximity
			move_agent!(particle, world)
		end
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
		:sticky_x => 0:1:20,
		:sticky_y => 0:1:20
	)

	plotkwargs = (
		agent_color=:green,
		agent_marker=(ag->wedge(ag,0.5)),
	)

	playground, = abmplayground( stabilisation; params, plotkwargs...)

	playground
end

end