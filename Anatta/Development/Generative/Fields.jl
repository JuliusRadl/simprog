#========================================================================================#
"""
	Fields

This module illustrates how to implement and handle the presence of Fields in an agent's world
that can react, diffuse and evaporate. Such flows provide the topological structure around
which biological form develops.

Author: Niall Palfreyman, March 2025.
"""
module Fields

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools, LinearAlgebra

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	SlimeMould

An extremely simple agent that just wanders around the world, endlessly depositing cyclic
Adenosine-MonoPhosphate (cAMP) into its environment.
"""
@agent struct SlimeMould(ContinuousAgent{2,Float64})
	preference::Float64		# how much to turn towards high/low cAMP concentration
	aversion::Bool			# if low cAMP is preferred
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	fields()

Create the Fields world model.
"""
function fields(;
	extent	= (50, 50),		# Dimensions of the world
	dt		= 0.1,			# Time step
	dcAMP	= 1.0,			# Deposition rate of cAMP
	ΛcAMP	= 0.0,			# Diffusion rate of cAMP
	αcAMP	= 0.0,			# Evaporation rate of cAMP
	preference = 0.001,
	aversion = false	
)
	# Create the Fields space, properties and model:
	model = StandardABM( SlimeMould, ContinuousSpace(extent; spacing=1.0);
        agent_step!, model_step!,
		properties = Dict(		# Model properties applying to all SlimeMoulds:
			:dt			=> dt,
			:dcAMP		=> dcAMP,
			:ΛcAMP		=> ΛcAMP,
			:αcAMP		=> αcAMP,
			:cAMP		=> zeros(Float64, extent),	# Heatmap of cAMP at each location in world
			:preference => preference,
			:aversion => aversion
		)
    )

	# Insert five SlimeMoulds at random points in the world, facing randomly:
	for _ in 1:5
		theta = 2π*rand()								# Random angle
		add_agent!( model, (cos(theta), sin(theta)), preference, aversion)
	end

	return model
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!(slimy, model)

The SlimeMould deposits cAMP at its current position.
"""
function agent_step!( slimy, model)
	idxHere = get_spatial_index( slimy.pos, model.cAMP, model)
	dt = model.dt

	model.cAMP[idxHere] += model.dcAMP * dt			# Deposit cAMP at current location

	# To-do: Implement niche-construction:
	#wiggle!(slimy, dt)								# Turn a little, then ...

	# velocity: a vector where its angle is the direction and its length (norm)
	# is its speed
	vel = slimy.vel
	speed = AgentTools.norm(vel)
	dir = vel ./ speed

	# gradient: the direction in which field value increases the most (also an
	# angle and a speed, but we are only interested in the angle)
	gradient = AgentTools.gradient(slimy.pos, model.cAMP, model)
	gradient = slimy.aversion ? gradient .* -1 : gradient
	gspeed = AgentTools.norm(gradient)
	# if gradient is 0 (no increase or decrease in any direction), keep original
	# direction to avoid dividing by zero
	gdir = gspeed > 0 ? gradient ./ gspeed : dir

	# create new direction as linear combination of direction and gradient
	# and adjust to original speed
	dirnew = LinearAlgebra.normalize((1 - slimy.preference) .* dir .+ slimy.preference .* gdir)
	slimy.vel = dirnew * speed
	# update agents facing angle
	move_agent!(slimy, model, 3)					# step forward with given speed
end

#-----------------------------------------------------------------------------------------
"""
	model_step!(model)

cAMP at every location diffuses through the world and evaporates by a constant fraction.
"""
function model_step!(model)
	dt = model.dt
	model.cAMP[:] = (
		(1-model.αcAMP*dt) .* diffuse4( model.cAMP, model.ΛcAMP*dt)
	)[:]
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Create a playground for visualising manipulation of cAMP field by a SlimeMould.
"""
function demo()
	params = Dict(												# Sliders:
		:ΛcAMP	=> 0:0.01:0.1,										# cAMP diffusion rate
		:αcAMP	=> 0:0.001:0.01,										# cAMP evaporation rate
		:preference => 0:0.001:1.0,
		:aversion => false:true
	)
	plotkwargs = (												# Plotting parameters:
		agent_color=:red, 											# SlimeMould's colour
		agent_marker=wedge,											# SlimeMould's marker
		#title = "Creation and manipulation of fields",				# Playground title
		heatarray=(model->model.cAMP), 								# Background map of cAMP
		heatkwargs = (colormap=:thermal, colorrange=(0.0,10.0)),	# cAMP intensity colours
        mdata=[(m->mean(m.cAMP))], mlabels=["Mean cAMP density"],	# Mean cAMP density data
	)

	playground, = abmplayground( fields; params, plotkwargs...)

	playground
end

end # ... end of module Fields
