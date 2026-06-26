#========================================================================================#
"""
	StructuralTuring

In this module, we replace the continuous flow dynamics of Turing systems by a discrete structural
algorithm that provides a reasonably good approximation of the original Turing dynamics. The
advantage of replacing dynamics by an algorithm is much higher computational performance, but we
should always be aware that such structural approximations might destroy important behaviours of
the original continuous system.

Author: Niall Palfreyman, March 2025.
"""
module StructuralTuring

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, Random, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	SlimeMould

A SlimeMould can in principle move, and secretes the chemical A(ctivator).
"""
@agent struct SlimeMould(ContinuousAgent{2,Float64})
	a_secretion_rate::Float64
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	structural_turing( kwargs)

Initialise the StructuralTuring algorithmic reaction-diffusion system.
"""
function structural_turing(;
	inner_radius_x = 2.5,
	inner_radius_y = 4.2,
	outer_radius_x = 6.0,
	outer_radius_y = 6.8,
	inhibition = 0.35,
	pPop = 0.01,
)
	extent = (101,101)
	properties = Dict(
		:inner_radius_x => inner_radius_x,
		:inner_radius_y => inner_radius_y,
		:outer_radius_x => outer_radius_x,
		:outer_radius_y => outer_radius_y,
		:inhibition => inhibition,
		:differentiation => zeros(Float64,extent...),
#		:differentiation => rand(Float64,extent...),
		:activators => fill(Int[],extent),
		:inhibitors => fill(Int[],extent),
		:pPop => pPop,
		:dt => 0.1,
	)

	st = StandardABM( SlimeMould, ContinuousSpace(extent,spacing=1.0);
		agent_step!, model_step!, properties
	)

	for i in 1:extent[1], j in 1:extent[2]
		st.activators[i,j] = nearby_patches( (i,j), st, (st.inner_radius_x,st.inner_radius_y))
		st.inhibitors[i,j] = setdiff(
			nearby_patches( (i,j), st, (st.outer_radius_x,st.outer_radius_y)),
			st.activators[i,j]
		)
	end

	if any(isempty.(st.inhibitors))
		println("Warning: outer_range must be greater than inner_range in both dimensions")
	end

	for _ in 1:pPop*prod(extent)
		add_agent!( SlimeMould, st, (0,0), 0.1)
	end

	st
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( slimy, st)

Express activator and inhibitor proteins.
"""
function agent_step!( slimy::SlimeMould, st::ABM)
	if slimy.a_secretion_rate > 0.0
		idx_here = get_spatial_index( slimy.pos, st.differentiation, st)
		st.differentiation[idx_here] += slimy.a_secretion_rate*st.dt
	end
end

#-----------------------------------------------------------------------------------------
"""
	model_step!(st)

Allow each location in the world to 
"""
function model_step!( st::ABM)
	for i in 1:length(st.differentiation)
		# Calculate total activatory influence on me ...
		activation = sum( st.differentiation[st.activators[i]]) -
			st.inhibition * sum( st.differentiation[st.inhibitors[i]])

		# ... then adjust differentiation level accordingly - either downwards towards zero:
		st.differentiation[i] *= (1-st.dt)
		if activation > 0
			# ... or else upwards:
			st.differentiation[i] += st.dt
		end
	end
end

#-----------------------------------------------------------------------------------------
"""
	nearby_patches( centre::Tuple, model::ABM, r::Tuple)

Return a Vector of all indices of model within range r from a centre location.
Note: This implementation is still quite inefficient and only works for 2d spaces.
"""
function nearby_patches( centre::Tuple, model::ABM, r::Tuple)
	# Set up neighbourhood Cartesian indices:
	rsq = r.^2
	rint = (s->round(Int,s)).(r)
	nbhd = [(i,j)
		for i in -rint[1]:rint[1], j in -rint[2]:rint[2]
			if (i^2/rsq[1] + j^2/rsq[2] ≤ 1.0)
	]
	setdiff!(nbhd,[(0,0)])						# We only want neighbours - not the centre

	# Displace (0,0)-based neighbourhood locations to the centre location:
	cint = (s->round(Int,s)).(centre)
	disc = map(nbhd) do pt
		cint .+ pt
	end

	# Extract linear matrix indices for the neighbourhood locations:
	ext = Int.(abmspace(model).extent)
#	ext = Int.(model.space.extent)
	ind = LinearIndices((w->1:w).(ext))

	map(disc) do tupl
		ind[mod(tupl[1]-1,ext[1])+1,mod(tupl[2]-1,ext[2])+1]
	end
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Create playground with two wave sources.
"""
function demo()
	params = Dict(
		:inner_radius_x => 0:0.1:10,
		:inner_radius_y => 0:0.1:10,
		:outer_radius_x => 0:0.1:10,
		:outer_radius_y => 0:0.1:10,
		:inhibition => 0:0.05:2,
		:pPop => 0:0.01:1,
	)
	plotkwargs = (
		agent_marker = :circle,
		agent_color = :red,
		agent_size = 15,
		heatarray = :differentiation,
		add_colorbar = false,
	)

	playground, = abmplayground( structural_turing; params, plotkwargs...)

	playground
end

end