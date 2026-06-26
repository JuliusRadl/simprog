#========================================================================================#
"""
	Boids

This module reproduces the work of Craig Reynolds in 1986 on flocking behaviour. It was later
used in cinema in the film Jurassic Park to generate the herd behaviour of dinosaurs running from
a predator.

Author: Niall Palfreyman, March 2025.
"""
module Boids

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Boid

The populating (flocking) agents in the Boids model.
"""
@agent struct Boid(ContinuousAgent{2,Float64})
	speed::Float64					# Speed of this Boid
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	boids( kwargs)

Create and initialise the Boids model.
"""
function boids(;
	n_boids = 50,					# Number of Boids in model
	speed = 1.0,					# Initial speed of Boids in model
	visual_range = 5.0,				# Visual eyesight range of Boids
	min_separation = 2.0,			# Preferred minimum separation between Boids
	cohesion_weight = 0.1,			# Priority of cohesion between Boids
	alignment_weight = 0.05,		# Priority of aligning velocity between Boids
	separation_weight = 0.25,		# Priority of separation between Boids
	extent = (100, 100),			# Extent of model space
)
	model = StandardABM( Boid, ContinuousSpace(extent; spacing=visual_range/2);
        agent_step!,
		properties = Dict(		# Model properties applying to all Boids:
			:n_boids			=> n_boids,
			:visual_range		=> visual_range,
			:min_separation		=> min_separation,
			:cohesion_weight	=> cohesion_weight,
			:alignment_weight	=> alignment_weight,
			:separation_weight	=> separation_weight,
		)
    )

	for _ in 1:n_boids
		theta = 2π*rand()								# Random angle
		add_agent!( model, (cos(theta), sin(theta)), speed)
	end

	return model
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( boid, flock)

This is the heart of the Boids model: It calculates precisely how each agent (Boid) behaves in
response to the behaviour of its nearest neighbours.
"""
function agent_step!(boid, flock)
	neighbour_ids = nearby_ids(boid, flock, flock.visual_range)
	n_nbrs = 0									# Number of visual neighbours
	align = separate = cohere = (0.0, 0.0)		# Upcoming adjustments to boid's behaviour

	for id in neighbour_ids
		# Apply Reynolds' flocking rules with respect to this neighbour:
		n_nbrs += 1
		nbrpos = flock[id].pos
		heading = get_direction(boid.pos,nbrpos,flock)	# Periodic-respecting displacement vector

		# SEPARATE to avoid colliding with neighbouring boid:
		if euclidean_distance(boid.pos, nbrpos, flock) < flock.min_separation
			separate = separate .- heading
		else
			# Align and cohere ONLY if neighbour is not too close!
			# ALIGN my flight to the trajectory of neighbouring boid:
			align = align .+ flock[id].vel
			# COHERE towards neighbouring boid:
			cohere = cohere .+ heading
		end
	end
	n_nbrs = max(n_nbrs, 1)

	# Calculate average velocity adjustments biased by priority:
	cohere		= cohere	./ n_nbrs .* flock.cohesion_weight
	align		= align		./ n_nbrs .* flock.alignment_weight
	separate	= separate	./ n_nbrs .* flock.separation_weight

	# Apply flocking rules to my current velocity:
	boid.vel = (boid.vel .+ cohere .+ align .+ separate) ./ 2
	boid.vel = boid.vel ./ norm(boid.vel)

	# Now move forward according to my newly adjusted velocity and speed:
	move_agent!(boid, flock, boid.speed)
end

#-----------------------------------------------------------------------------------------
function demo()
	params = Dict(
		:min_separation		=> 0:0.5:10,
		:cohesion_weight	=> 0:0.05:1,
		:alignment_weight	=> 0:0.05:1,
		:separation_weight	=> 0:0.05:1,
	)

    plotkwargs = (
        # To-do: Specify plotting keyword arguments
        agent_marker    = wedge,
        agent_color     = multicoloured,
    )

	playground, = abmplayground( boids; params, plotkwargs...)

	playground
end

end