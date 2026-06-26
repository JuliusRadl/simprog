#========================================================================================#
"""
	Replicators

Module Replicators: A model of an exponentially replicating population.

Author: Niall Palfreyman, 04/09/2022
"""
module Replicators

# Externally callable methods of Replicators:
#export Replicator

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
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
"""
	demo()

Demonstrate use of the Replicators module.
"""
function demo()
	println("\n============ Demonstrate Replicators: ===============")
	println("An exponential population of replicators from t=0-5 generations:")
	repl = Replicator( [0,1,2,3,4,5], 1, [0,0,0,0,0,0])
	display( repl)
	println()

	println("Run population with initial size x0=1 and growth constant mu=1:")
	# run!(repl,1)
	display( repl)
	println()

	println("Run population with initial size x0=1 and growth constant mu=2:")
	# display( run!(repl,1,2))
	println()

	println("Run population with initial size x0=3 and growth constant mu=1:")
	# display( run!(repl,3))
	println()
end

end		# ... of module Replicators

