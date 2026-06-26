#========================================================================================#
"""
	PSO

This model demonstrates how a Particle Swarm can solve a minimisation problem - and also the major
difficulty of minimisation: suboptimisation.
	
Author: Niall Palfreyman, March 2025.
"""
module PSO

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Ant

Ants run around the world seeking to minimise the value of some objective function u, but unlike
the Scouts of the Suboptimisation model, each Ant has a memory of their best-yet objective value.
"""
@agent struct Ant( ContinuousAgent{2,Float64})
	speed::Float64
	memory::Float64					# Current lowest value of u yet found in my travels
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	pso()

Initialise a PSO model.
"""
function pso(;
	difficult=false,
	temperature=0.002,
)
	extent = (99,99)
	pPop = 0.01			# Proportion of populated locations

	properties = Dict(
		:objective => (difficult ? dejong2(extent) : valleys(extent)),
		:pPop => pPop,
		:difficult => difficult,
		:temperature => temperature,
		:tolerance => 0.1,
		:mass_centre => 0.5 .* extent
	)

	world = StandardABM( Ant, ContinuousSpace( extent; spacing=1.0);
		model_step!, agent_step!, properties
	)
	
	for _ in 1:pPop*prod(extent)
		# Baby Ant has unit speed and ridiculously high memory of objective function:
		θ = 2π * rand()
		add_agent!( world, (cos(θ),sin(θ)), 1.0, 1e301)
	end

	world
end

#-----------------------------------------------------------------------------------------
"""
	model_step!( world)

Recalculate mass_centre of all Ants in the world.
"""
function model_step!(model)
	mc = (0.0,0.0)
	for ant in allagents(model)
		mc = mc .+ ant.pos
	end
	model.mass_centre = mc ./ nagents(model)

	return model
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( ant, world)

This is the stabilising procedure for Ants: Wiggle randomly if objective value u is rising, but
maintain direction of travel if u is decreasing. Also, slow down if other agents are around.
Finally, add some random motion - especially if you remember you have previously encountered a
lower objective value than the current one.
"""
function agent_step!( ant, world)
	nbr_range = 1.1
	u = world.objective
	prev_u = u[get_spatial_index( ant.pos, u, world)]

	if rand() < ant.speed
		move_agent!( ant, world, nbr_range)
	end
	curr_u = u[get_spatial_index( ant.pos, u, world)]
	nbrs = collect(nearby_agents( ant, world, nbr_range))

	# Remember lower values of u - otherwise spread dissatisfaction:
	if curr_u < ant.memory + world.tolerance
		ant.memory = curr_u
		accelerate!(ant,false)
	else
		if length(nbrs) != 0
			accelerate!(rand(nbrs))
		end
	end

	# Wiggle if objective is increasing:
	if curr_u > prev_u
		wiggle!(ant,pi)
	end

	# Be sociable - slow down for neighbours:
	if length(nbrs) != 0
		accelerate!(ant,false)
	end

	# Add in some random motion:
	if rand() < world.temperature
		accelerate!(ant)
	end
end

#-----------------------------------------------------------------------------------------
"""
	accelerate!( ant::Ant, positively::Bool=true)

If ant is positively accelerating, increase speed; otherwise slow down.
"""
function accelerate!( ant::Ant, positively::Bool=true)
	inertia = 0.5
	if positively
		ant.speed = inertia*ant.speed + (1-inertia)
	else
		ant.speed *= inertia
	end
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Demonstrate the PSO model.
"""
function demo()
	params = Dict(
		:difficult => false:true,
		:temperature => 0.0:0.002:1.0
	)

	plotkwargs = (
		agent_colour=:lime,
		agent_size=15,
		add_colorbar=false,
		heatarray=:objective,
	)

	playground,obs = abmplayground( pso; params, plotkwargs...)

	# Whenever mass_centre is changed, update its representation on the heatmap:
	mcx = lift( (wld->wld.mass_centre[1]), obs.model)
	mcy = lift( (wld->wld.mass_centre[2]), obs.model)
	scatter!( mcx, mcy, color=:darkgreen, markersize=30)

	playground
end

end