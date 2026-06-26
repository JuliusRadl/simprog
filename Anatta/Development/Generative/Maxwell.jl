#========================================================================================#
"""
	Maxwell

This module encapsulates the simulation of an electromagnetic field defined by Maxwell's equations.
It is capable of simulating the propagation of electromagnetic waves in a 2D space.

Author: Niall Palfreyman, March 2025
"""
module Maxwell

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Source

An electrical field source (possibly including a time-varying component).
"""
@agent struct Source(ContinuousAgent{2,Float64})
	q0::Float64			# Initial charge
	f::Float64			# Frequency of oscillation
	phi::Float64		# Phase of oscillation
	q::Float64			# Current charge of the source
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	maxwell( kwargs)

Initialise Maxwell model initially empty of Sources.
"""
function maxwell(;
	d			= 10.0,			# Spacing between sources
	attenuation = 0.4,			# Attenuation factor for wave propagation (0.0 = no attenuation)
)
	extent = (100,100)			# Extent of model space
	world = ContinuousSpace(extent, spacing=1.0)
	properties = Dict(
		:E => zeros(Float64,extent),
		:B => zeros(Float64,extent),
		:c => 1.0,
		:d => d,
		:t => 0.0,
		:dt => 0.004,
		:dxy => 0.1,
		:attenuation => attenuation,
	)

	maxi = StandardABM( Source, world; agent_step!, model_step!, properties)
	halfway = minimum(extent)÷2
	addsource!( maxi, (halfway,halfway-d/2))
	addsource!( maxi, (halfway,halfway+d/2))

	maxi
end

#-----------------------------------------------------------------------------------------
"""
	addsource!( maxi, pos, amplitude, frequency, phase)

Add a Source to the Maxwell model with the given specification, initialised to time zero.
"""
function addsource!( maxi, pos; vel=(0.0,0.0), q0=1.0, f=1.0, phi=0.0)
	add_agent!( pos, maxi, vel, q0, f, phi, q0*sin(phi))
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( src, maxi)

Source agent updates local E-field and moves to new position. Local E-field updates are then
propagated outwards by the wave equation implemented in model_step!.
"""
function agent_step!( src::Source, maxi::ABM)
	maxi.E[get_spatial_index( src.pos, maxi.E, maxi)] = src.q =
		src.q0 * sin(src.phi + 2pi*src.f*maxi.t)
	move_agent!( src, maxi, maxi.dt)
end

#-----------------------------------------------------------------------------------------
"""
	model_step!(maxi)

Approximate a wave solution of Maxwell's equations in 2D space by updating the B-field according
to the Laplacian of the E-field, and then updating the E-field according to the B-field, where
the Laplacian operator is approximated by the 2-D curvature of the E-field:
	d2E/dt2 = c^2 * del^2(E)
"""
function model_step!( maxi::ABM)
	# First calculate first-order Euler change in B due to current E configuration ...
	maxi.B += maxi.dt * maxi.c^2 * (
		0.25(
			circshift(maxi.E,( 1,0)) + circshift(maxi.E,(0, 1)) +
			circshift(maxi.E,(-1,0)) + circshift(maxi.E,(0,-1))
		) - maxi.E
	) / maxi.dxy^2
	maxi.B *= (1-maxi.attenuation*maxi.dt)			# Attenuate the wave

	# ... then calculate second-order change in E due to current B configuration:
	maxi.E += maxi.dt * maxi.B
	maxi.t += maxi.dt
end

#-----------------------------------------------------------------------------------------
"""
	charge( src::Source)

Adapt colour of the source agent according to its charge.
"""
charge( src::Source) = (src.q >= 0 ? :yellow : :black)

#-----------------------------------------------------------------------------------------
"""
	demo()

Create playground with two wave sources.
"""
function demo()
	params = Dict(
		:attenuation	=> 0:0.1:5,
		:d				=> 4:20.0,
	)
	plotkwargs = (
		agent_color = charge,
		agent_size = 20,
		heatarray = :E,
		add_colorbar = false,
		heatkwargs = (
			colorrange = (-1,1),
			colormap = cgrad(:oslo),
		)
	)

	playground, = abmplayground( maxwell; params, plotkwargs...)

	playground
end

end