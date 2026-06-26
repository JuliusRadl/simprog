#========================================================================================#
"""
	Interactors

Module Interactors: A model of evolutionary game theory.

Author: Julius Radl, 2025-12-08
"""
module Interactors

include( "Simplex.jl")

using GLMakie, LinearAlgebra

export Interactor, simulate!, plot3!

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Interactor

A Interactor represents the time-evolution of population frequencies
considering interactive fitness
"""
struct Interactor
	resolution::Int				# Resolution of timescale
	t::Vector{Float64}			# Timescale
	x::Vector{Vector{Float64}}	# Time-series of three population types
	payoffs::Matrix{Float64}	# 3x3 Matrix describing encounter payoffs


 	function Interactor(payoffs, r::Vector{Float64} = [1., 1., 1.])
		ngenerations = 1000
		new(ngenerations, zeros(Float64,ngenerations+1),
		Vector{Vector{Float64}}(undef,ngenerations+1), payoffs)
	end
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	plot3!( axis, interactor)

Display the simulation results in interactor as a simplex.
"""
function plot3!( axis::Axis, int::Interactor)
	Simplex.plot3!( axis)			# Plot the axes
	Simplex.plot3!( axis, int.x)	# Plot the trajectory
end

"""
	simulate!( interactor, x0, T)

Simulate the interactor dynamics for T ticks, starting from initial state x0.
"""
function simulate!( int::Interactor, x0::Vector{Float64}, T::Real)
	dt = T/int.resolution;			# Full time-step
	dt2 = dt/2;						# RK2 half time-step
	
	# Set up time scale and initial, normalized value of population:
	int.t[:] = 0:dt:T
	int.x[1] = x0 / sum(x0)
	
	# Calculate population trajectory:
	for step = 1:(int.resolution)
		# Legibility
		x = int.x[step]
		A = int.payoffs

		# Perform RK2 half-step:
		R = x' * A * x
		dxdt = x .* (A * x .- R)
		half_step = x + x .* dxdt * dt2
		
		# Perform RK2 full-step:
		new_R = half_step' * A * half_step
		new_dxdt = half_step .* (A * half_step .- new_R)
		int.x[step+1] = x + x .* new_dxdt * dt

		# Normalize frequencies
		int.x[step+1] = int.x[step+1] / sum(int.x[step+1])
		end
end

"""
	unittest()

Unit-test the Interactors module.
"""
function unittest()
	println("\n============ Unit test Interactors: ===============")

	# Set up plot:
	fig = Figure()
	ax1 = Axis(fig[1, 1])
	ax2 = Axis(fig[1, 2])
	ax3 = Axis(fig[2, 1])
	ax4 = Axis(fig[2, 2])
	
	# Cyclic Dominance: Jeder darf mal viele sein
	payoffs = [0 1 -1; -1 0 1; 1 -1 0]
	int = Interactor(payoffs)
	simulate!( int, [10.,2.,1.], 1000)
	plot3!(ax1,int).parent

	# Cyclic Dominance in klein
	payoffs2 = [0 1 -1; -1 0 1; 1 -1 0]
	int2 = Interactor(payoffs2)
	simulate!( int2, [2., 1., 1.], 1000)
	plot3!(ax2,int2)

	# Snakes, Doves, Hawks
	payoffs3 = [1 4 0; 0 5 0; 0 4 1]
	int3 = Interactor(payoffs3)
	simulate!( int3, [1., 3., 1.], 1000)
	plot3!(ax3, int3).parent
end

end		# ... of module Interactors