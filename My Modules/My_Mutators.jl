#========================================================================================#
"""
	Mutators

Module Mutators: A Model of mutations and population frequency development 
in three quasi-species

Author: Julius Radl, 2025-12-07
"""
module Mutators

include( "Simplex.jl")

using GLMakie, LinearAlgebra

export Mutator, simulate!, plot3!

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Mutator

A Mutator represents the time-evolution of a mutation model
"""
struct Mutator
	r::Vector{Float64}				# Vector of three type fitnesses
	resolution::Int				# Resolution of timescale
	t::Vector{Float64}			# Timescale
	x::Vector{Vector{Float64}}	# Time-series of three population types
	mut_mat::Matrix{Float64}	# 3x3 stochastic matrix describing mutation rates


 	function Mutator(mut_mat, r::Vector{Float64} = [1., 1., 1.])
		ngenerations = 1000
		# TO DO: Normierung streichen? Im Skript ist R = r = 1, nicht 0.33
		new(r/sum(r), ngenerations, zeros(Float64,ngenerations+1),
		Vector{Vector{Float64}}(undef,ngenerations+1), mut_mat)
	end
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	plot3!( axis, mutator)

Display the simulation results in mutator as a simplex.
"""
function plot3!( axis::Axis, mut::Mutator)
	Simplex.plot3!( axis)			# Plot the axes
	Simplex.plot3!( axis, mut.x)	# Plot the trajectory
end

"""
	simulate!( mutator, x0, T)

Simulate the mutator dynamics for T ticks, starting from initial state x0.
"""
function simulate!( mut::Mutator, x0::Vector{Float64}, T::Real)
	dt = T/mut.resolution;			# Full time-step
	dt2 = dt/2;						# RK2 half time-step
	
	# Set up time scale and initial, normalized value of population:
	mut.t[:] = 0:dt:T
	mut.x[1] = x0 / sum(x0)
	
	# Calculate population trajectory:
	for step = 1:(mut.resolution)
		# Legibility
		r = mut.r
		x = mut.x[step]
		Q = mut.mut_mat

		# Perform RK2 half-step:
		dxdt = Q * (x .* r) - (x' * r) * x
		half_step = x + x .* dxdt * dt2
		
		# Perform RK2 full-step:
		new_dxdt = Q * (half_step .* r) - (half_step' * r) * half_step
		mut.x[step+1] = x + x .* new_dxdt * dt

		# Normalize frequencies
		mut.x[step+1] = mut.x[step+1] / sum(mut.x[step+1])
		end
end

"""
	unittest()

Unit-test the Mutators module.
"""
function unittest()
	println("\n============ Unit test Mutators: ===============")

	# Set up plot:
	fig = Figure()
	ax1 = Axis(fig[1, 1])
	ax2 = Axis(fig[1, 2])
	ax3 = Axis(fig[2, 1])
	ax4 = Axis(fig[2, 2])
	
	# Cyclic Mutation with even fitness values
	mut_mat = [[0.3, 0.3, 0.4] [0.4, 0.3, 0.3] [0.3, 0.4, 0.3]]
	mut = Mutator(mut_mat)
	simulate!( mut, [10.,2.,1.], 1000)
	plot3!(ax1,mut)

	# Cyclic Mutation, Type 3 most fitness
	mut_mat2 = [[.9, .1, .0] [.0, .9, .1] [.1, .0, .9]]
	mut2 = Mutator(mut_mat2, [.1, .2, .7])
	simulate!( mut2, [10.,2.,1.], 100)
	plot3!(ax2,mut2)

	# Type 1 most stable but least fitness
	mut_mat3 = [[0.9, .1, 0] [.3, .6, .1] [.4, .1, .5]]
	mut3 = Mutator(mut_mat3, [.1, .2, .7])
	simulate!( mut3, [1., 2., 10.], 100)
	plot3!(ax3,mut3)

	# No Mutation, No Selection, No Change
	mut_mat4 = [[1, 0, 0] [0, 1, 0] [0, 0, 1]]
	mut4 = Mutator(mut_mat4)
	simulate!( mut4, [10.,2.,1.], 100)
	plot3!(ax4,mut4).parent
end

end		# ... of module Mutators