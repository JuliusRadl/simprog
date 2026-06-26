#========================================================================================#
"""
	Ant search algorithms

A simulation to investigate manipulation of the environment to enable problem-solving.
	
Author: Niall Palfreyman, March 2025.
"""
module AntSearch

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Ant

Ants forage for food, which they bring back to the nest. Their state is defined by how much
they are carrying.
"""
@agent struct Ant(ContinuousAgent{2,Float64})
	carrying::Float64										# How much food am I carrying?
end

#-----------------------------------------------------------------------------------------
"""
	FoodSource

FoodSources do nothing in this model other than provide food at a certain constant capacity.
"""
@agent struct FoodSource(ContinuousAgent{2,Float64})
	capacity::Float64
end

#-----------------------------------------------------------------------------------------
"""
	Nest

Nests are pretty boring: they just sit there as a receptacle for foraged food.
"""
@agent struct Nest(ContinuousAgent{2,Float64})
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	ant_search()

Initialise the AntSearch simulation.
"""
function ant_search(;
	diff_rate = 0.6,
	evap_rate = 0.05,
	upper_right_capacity = 2.0,
	lwb_tolerance = 0.05,
	upb_tolerance = 1.0,
)
	extent = (70,70)
	pPop = 0.02
	population = round(Int,pPop*prod(extent))
	nest_x = extent[1]/2

	# Calculate nest_pheromone concentration centred on nest_pos:
	y_distance = (1:minimum(extent)) .- nest_x
	squ_y_distance = repeat(y_distance.^2, 1, minimum(extent))
	nest_distance = sqrt.(squ_y_distance + squ_y_distance')
	max_distance = findmax(nest_distance)[1]
	nest_phero = max_distance .- nest_distance

	properties = Dict(
		:carry_phero => zeros(Float64,extent),			# Strength of carrying pheromone here
		:nest_phero => nest_phero,						# Strength of nest pheromone here
		:sniff_radius => 3.0,							# Effective radius of nest and food sources
		:diff_rate => diff_rate,						# Diffusion rate of carry_phero
		:evap_rate => evap_rate,						# Evaporation rate of carry_phero
		:secrete_rate => 10.0,							# At what rate do ants secrete carry_pheno?
		:upper_right_capacity => upper_right_capacity,	# Capacity of food source 3
		:nest_pos => (nest_x,nest_x),					# Position of nest
		:food_source => Vector{FoodSource}(undef,3),	# The three food sources
		:lwb_tolerance => lwb_tolerance,				# Ant can't detect carry_phero below this
		:upb_tolerance => upb_tolerance,				# Ant ignores carry_phero above this
	)

	asearch = StandardABM(
		Union{Ant,FoodSource,Nest},
		ContinuousSpace(extent,spacing=1.0);
		properties, warn=false,
		model_step!, agent_step!
	)

	# Create the nest:
	add_agent!( asearch.nest_pos, Nest, asearch, (0,0))
	
	# Create food sources:
	asearch.food_source[1] = add_agent!( (0.1extent[1],0.7extent[2]), FoodSource, asearch, (0,0), 1.0)
	asearch.food_source[2] = add_agent!( (0.6extent[1],0.1extent[2]), FoodSource, asearch, (0,0), 2.0)
	asearch.food_source[3] = add_agent!( 0.9.*extent, FoodSource, asearch, (0,0), upper_right_capacity)

	# Scatter the Ants in search mode (i.e., carrying nothing):
	map(2π*rand(population)) do θ
		add_agent!( Ant, asearch, (cos(θ),sin(θ)), 0.0)
	end

	asearch
end

#-----------------------------------------------------------------------------------------
"""
	model_step!( asearch)

Diffuse and evaporate the carrying_pheromone of the Ants.
"""
function model_step!( asearch::ABM)
	asearch.carry_phero[:] = (
		(1-asearch.evap_rate) .* diffuse4( asearch.carry_phero, asearch.diff_rate)
	)[:]
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( ant, asearch)

One step in the life of an Ant engaged in searching for food.
"""
function agent_step!( ant::Ant, asearch)
	if ant.carrying > 0
		go_to_nest!(ant,asearch)
	else
		look_for_food!(ant,asearch)
	end

	turn!(ant,rand()-0.5)						# Wiggle ...
	move_agent!(ant,asearch,1.0)				# ... and step forward
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( fsnest, asearch)

FoodSources and Nests do nothing.
"""
function agent_step!( fsnest::Union{FoodSource,Nest}, asearch)
	dummystep(fsnest,asearch)
end

#-----------------------------------------------------------------------------------------
"""
	go_to_nest!(ant,asearch)

Carry your food to the nest and dump it.
"""
function go_to_nest!( ant::Ant, asearch)
	if euclidean_distance(ant.pos,asearch.nest_pos,asearch) < asearch.sniff_radius
		# We've arrived!
		ant.carrying = 0.0
		turn!(ant,pi)
	else
		# Still looking for nest ...
		idx_here = get_spatial_index( ant.pos, asearch.nest_phero, asearch)
		idx_fwd  = get_spatial_index( normalize_position(ant.pos.+ant.vel,asearch),
			asearch.nest_phero, asearch
		)
		if asearch.nest_phero[idx_fwd] < asearch.nest_phero[idx_here]
			turn!(ant,pi)
		end
		asearch.carry_phero[idx_here] += asearch.secrete_rate * ant.carrying
	end
end

#-----------------------------------------------------------------------------------------
"""
	look_for_food!(ant,asearch)

Wander around looking for food. If you find some, pick it up.
"""
function look_for_food!( ant::Ant, asearch)
	foodsrc = nothing
	for fs in asearch.food_source
		if euclidean_distance(ant.pos,fs.pos,asearch) < asearch.sniff_radius
			foodsrc = fs
			break
		end
	end

	if foodsrc !== nothing
		# We've found food:
		ant.carrying = foodsrc.capacity
		wiggle!(ant)
	else
		# Still looking for food ...
		idx_here = get_spatial_index( ant.pos, asearch.carry_phero, asearch)

		if asearch.lwb_tolerance < asearch.carry_phero[idx_here] < asearch.upb_tolerance
			# Carrying pheromone is distinctive, but not overwhelming:
			idx_fwd = get_spatial_index( normalize_position(ant.pos.+ant.vel,asearch),
				asearch.carry_phero, asearch
			)
			if asearch.carry_phero[idx_fwd] < asearch.carry_phero[idx_here]
				wiggle!(ant)
			end
		end
	end
end

#-----------------------------------------------------------------------------------------
"""
Assorted visual functionlets.
"""
ashape(ant::Ant) = :utriangle
ashape(nest::AbstractAgent) = :circle
acolour(ant::Ant) = ant.carrying>0 ? :orange : :red
acolour(nest::Nest) = :black
acolour(fs::FoodSource) = :lime
asize(ant::Ant) = 20
asize(nest::Nest) = 50
asize(fs::FoodSource) = 50

#-----------------------------------------------------------------------------------------
"""
	demo()

Demonstrate an ant search algorithm.
"""
function demo()
	params = Dict(
		:diff_rate => 0:0.01:1,
		:evap_rate => 0:0.01:0.2,
		:upper_right_capacity => 2:0.1:4,
		:lwb_tolerance => 0:0.01:0.1,
		:upb_tolerance => 0:0.1:5,
	)
	plotkwargs = (
		agent_marker = ashape,
		agent_color = acolour,
		agent_size = asize,
		heatarray = :carry_phero,
		add_colorbar = false,
	)

	playground, = abmplayground( ant_search; params, plotkwargs...)

	playground
end

end