#========================================================================================#
"""
	Replicators

Module Replicators: A model of an exponentially replicating population.

Author: Niall Palfreyman, 04/09/2022
"""
module Replicators

# Externally callable methods of Replicators:
export Replicator, run!

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Replicator

Replicator is a concrete data-type. That is, it is not just an abstract name that guides multiple
dispatch, but contains concrete data that it stores in memory. It encapsulates (i.e., hides from
the outside world) three concrete items of data: a time-scale t in time-steps dt, and a
corresponding time-series x representing the size of the replicator population over this
time-scale. The initial time-series consists solely of zeros.
"""
struct Replicator
	t::Vector{Real}			# The simulation time-scale
	dt::Real				# The simulation time-step
	x::Vector{Real}			# The population time-series

	function Replicator(duration::Real, dt=1)
		timescale = 0.0:dt:duration
		new(timescale, dt, zeros(Float64, length(timescale)))
	end
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------

"""
	run!( replicator, x0, mu=1.0)

Simulate the exponential growth of a Replicator population, starting from the initial
value x0, and with specific growth constant mu.
"""
function run!(repl::Replicator, x0::Real, mu::Real=1.0)
	repl.x[1] = x0                  # Set initial value of the population

	# Euler's method

	# for i in 2:length(repl.t)
	# 	# Perform Euler step:
	# 	repl.x[i] = repl.x[i-1] + repl.dt * mu * repl.x[i-1]
	# end

	# Runge-Kutta-2

	dt_2 = repl.dt/2                                # dt2 is one half-timestep
    for i in 2:length(repl.t)
        # Perform Runge-Kutta-2 step:
        x_2 = repl.x[i-1] + mu*dt_2*repl.x[i-1]     # Calculate new x halfway thru step
        repl.x[i] = repl.x[i-1] + mu*repl.dt*x_2    # Use x2 as better approximation
    end

	repl                                                                                    # Return the Replicator
end

"""
	demo()

Demonstrate use of the Replicators module.
"""
function demo()
	println("\n============ Demonstrate Replicators: ===============")
	println("An exponential population of replicators from t=0-5 generations:")
	repl = Replicator(5,1)
	display( repl)
	println()

	println("Run population with initial size x0=1 and growth constant mu=1:")
	run!(repl,1)
	display( repl)
	println()

	println("Run population with initial size x0=1 and growth constant mu=2:")
	display( run!(repl,1,2))
	println()

	println("Run population with initial size x0=3 and growth constant mu=1:")
	display( run!(repl,3))
	println()
end

end		# ... of module Replicators

