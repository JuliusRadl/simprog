#========================================================================================#
"""
	Cooperators

Module Cooperators: This model simulates the evolution of cooperation between rats. In
particular, it implements a population frequency model containing n rat types, each
characterised by its own PD strategy.

Author: Julius Radl, 2025-12-08
"""
module Cooperation

include( "Simplex.jl")

using GLMakie, LinearAlgebra, Printf

export Rats, simulate!, plot3!

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Rats()

A population of different rat types with different strategies.
"""
mutable struct Rats
	resolution::Int				# Resolution of timescale
	t::Vector{Float64}			# Timescale
	x::Vector{Vector{Float64}}	# Time-series of type frequencies
	S::Matrix{Float64}			# Strategy matrix
	A_pd::Matrix{Float64}		# Prisoners' dilemma payoff matrix
	A::Matrix{Float64}			# Long term payoff matrix

	function Rats(nStrategies = 5, A_pd = [4. 0.; 5. 1.])
		resolution = 1000
		t = Vector{Float64}(undef, resolution+1)
		x = Vector{Vector{Float64}}(undef, resolution+1)
		S = rand(Float64, nStrategies, 2)					# Randomly generate strategy matrix
		x0_raw = rand(Float64, nStrategies)					# Randomize initial type frequencies
		x[1] = x0_raw ./ sum(x0_raw)
		A = calc_longterm_payoff(A_pd, S)
		
		new(resolution, t, x, S, A_pd, A)
	end
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	calc_longterm_payoff(A_pd)

Helper function for calculating the long term payoff matrix for a given
Prisoner's Dilemma Payoff Matrix (2 x 2) and a vector of strategies (n x 2).
"""
function calc_longterm_payoff(A_pd::Matrix{Float64}, S::Matrix{Float64})
	
	Te = A_pd[2, 1]	# Temptation to defect
	Re = A_pd[1, 1]	# Reward for mutual cooperation
	Pu = A_pd[2, 2]	# Punishment for mutual defection
	Su = A_pd[1, 2]	# Sucker's payoff

	p = S[:, 1]		# Reciprocation
	q = S[:, 2]
	r = p .- q

	s_ij = (r .* q' .+ q) ./ (1 .- (r * r'))	# Prob. that S[i] cooperates with S[j] long-term
	s_ji = s_ij'								# The reverse

	
	Te_part = Te * ((1 .- s_ij) .* s_ji)
	Re_part = Re * (s_ij .* s_ji)
	Pu_part = Pu * ((1 .- s_ij) .* (1 .- s_ji))
	Su_part = Su * (s_ij .* (1 .- s_ji))

	A = Te_part + Re_part + Pu_part + Su_part	# Adding matrices of the same size is no problem

end

#-----------------------------------------------------------------------------------------
"""
	plot3!( axis, cooperator)

Display the simulation results in cooperator as a simplex.
"""
function plot3!( axis::Axis, rats::Rats)
	Simplex.plot3!( axis)					# Plot the axes
	Simplex.plot3!( axis, rats.x[:][1:3])	# Plot the trajectory
											# of the first three strategies
end

#-----------------------------------------------------------------------------------------
"""
	simulate!(rats, nGenerations)

Iterate the replicator equation on the rat population over T ticks.
"""
function simulate!( rats::Rats, T::Real)
	dt = T/rats.resolution;			# Full time-step
	dt2 = dt/2;						# RK2 half time-step
	
	# Set up time scale
	rats.t[:] = 0:dt:T
	
	# Calculate population trajectory:
	for step = 1:(rats.resolution)
		# Legibility
		x = rats.x[step]
		A = rats.A

		# Perform RK2 half-step:
		R = x' * A * x
		dxdt = x .* (A * x .- R)
		half_step = x + x .* dxdt * dt2
		
		# Perform RK2 full-step:
		new_R = half_step' * A * half_step
		new_dxdt = half_step .* (A * half_step .- new_R)
		rats.x[step+1] = x + x .* new_dxdt * dt

		# Normalize frequencies
		rats.x[step+1] = rats.x[step+1] / sum(rats.x[step+1])
		end
end

#-----------------------------------------------------------------------------------------
"""
	mutate!(rats)

Replace the worst performing strategy with a new, random strategy.
"""
function mutate!(rats)

	# Legibility
	x = rats.x[end]

	# Find number of Strategies
	nStrategies = length(x)
	tenPercent = ceil(Integer, nStrategies/10)

	# Find worst strategies
	order = sortperm(x)
	toBeEliminated = order[1:tenPercent]
	# Generate new strategies
	rats.S[toBeEliminated, :] = rand(Float64, tenPercent, 2)

	# Recalculate Payoff matrix
	rats.A = calc_longterm_payoff(rats.A_pd, rats.S)
	
end

#-----------------------------------------------------------------------------------------
"""
	cooperation_axis(parent, xVector, yVector)

Plots two vectors of equal length as x- and y-components of points.
"""
function cooperation_axis(parent, pVector, qVector)
	# Create new axis
	ax = Axis(
		parent,
		xlabel = "Reciprocation",
		ylabel = "Forgiveness",
		title = "Trajectory of Reciprocation against Forgiveness over time",
		limits = (0., 1., 0., 1.),
		xticks = 0:0.1:1,
		yticks = 0:0.1:1
	)

	# Create the plot
	scatterlines!(ax, pVector, qVector, color = :blue, linewidth = 2, markercolor = :red, markersize = 10)
end

#-----------------------------------------------------------------------------------------
"""
	displayStrategies(rats)

Display the top n strategies in a legible manner.
"""
function displayStrategies(rats::Rats, nTopStrategies::Int64)
	order = sortperm(rats.x[end], rev=true)[1:nTopStrategies]

	println("Platz:\tp:\tq:\tx:")
	for i in 1:5
		j = order[i]
		@printf("%d\t%.2f\t%.2f\t%.2f\n", i, rats.S[j, 1], rats.S[j, 2], rats.x[end][j])
	end
	println()
end

#-----------------------------------------------------------------------------------------
"""
	unittest()

Unit-test the Cooperation module for development purposes.
"""
function unittest()
	println("\n============ Unit test Cooperation: ===============")

end

#-----------------------------------------------------------------------------------------
"""
	demo()

Demonstrate how to use the Cooperation module.
"""

function demo()
	nStrategies = 5
	topStrategies = 5
	nGenerations = 20
	nMutations = 1000
	rats = Rats( nStrategies)           # Create n-strategy population

	averages = Matrix{Float64}(undef, nMutations, 2)

	for mut in 1:nMutations
		simulate!( rats, nGenerations)					# Iterate replicator equation for nGenerations
		x = rats.x[end]
		p = rats.S[:, 1]
		q = rats.S[:, 2]
		averages[mut, :] = [dot(x, p), dot(x, q)]
		displayStrategies( rats, topStrategies)			# Display best strategies
		mutate!(rats)									# Replace worst strategy with new, random strategy
	end;

	simulate!( rats, nGenerations)
	displayStrategies( rats, topStrategies)
	fig = Figure()
	cooperation_axis(fig[1, 1], averages[:, 1], averages[:, 2])
	fig
end

end		# ... of module Cooperators