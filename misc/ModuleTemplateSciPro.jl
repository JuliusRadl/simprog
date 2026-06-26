#========================================================================================#
"""
	FunGraphics

This is a module designed to introduce you to the fun of using graphics in julia. It includes
two main elements: the user-defined type HillTF (that is, a Hill-function transcription factor),
and a fun_graphics function which gives you a chance to think about how graphics commands are
put together. Have fun! :)

Author: Niall Palfreyman, 26/05/2023
"""
module FunGraphics

using Observables, GLMakie

import Base.==
export animate_hill, expression, HillTF, fun_graphics

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	HillTF

This data type encapsulates the chemically important characteristics of a range of concentrations
of a Hill transcription factor, a half-response constant K and a cooperation level n.
"""
struct HillTF
	range::Vector{Float64}
	K::Float64
	n::Float64

	function HillTF( range::AbstractVector, K=2.0, n=1.0)
		@assert (K > 0.0) "K must be greater than zero"
		new(collect(range),K,n)
	end
end

#-----------------------------------------------------------------------------------------
# Module data:
#-----------------------------------------------------------------------------------------
"""
	scrambled_code

The scrambled code used for fun_graphics.
"""
const scrambled_code = [
	"    markersize=300*abs.(colours),"
	"fig, ax, plt = scatter( xdata, ydata;"
	"    figure=(; resolution=(600,400))"
	"limits!(0.0,1.0,0.0,1.0);"
	"    color=colours, label=\"Bubbles\", colormap=:plasma,"
	"Colorbar( fig[1,2], plt, height=Relative(3/4));"
	"fig"
	");"
	"xdata = rand(50); ydata = rand(50); colours = rand(50);"
	"Legend( fig[1,2], ax, valign=:top);"
	"    axis=(; aspect=DataAspect()),"
]

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	==( tf1::HillTF, tf2::HillTF) :: Bool

Equality between two HillTFs.
"""
function ==( tf1::HillTF, tf2::HillTF)
	tf1.K==tf2.K && tf1.n==tf2.n && tf1.range==tf2.range
end

"""
	hill( tf::Real, K, n)

Return the Hill chemical reaction rate corresponding to a transcription factor concentration tf,
a half-response constant K and a cooperation level n.
"""
function hill( tf::Real, K=2.0, n=1.0)
	abs_n = abs(n)
	if abs_n != 1
		K = K^abs_n
		tf = tf^abs_n
	end

	(n >= 0) ? (tf/(K+tf)) : (K/(K+tf))
end

#-----------------------------------------------------------------------------------------
"""
	expression( substrate::HillTF)

Return the protein expression rates corresponding to a range of concentrations of a transcription
factor HillTF, a half-response constant K and a cooperation level n.
"""
function expression( tf::HillTF)
	map(tf.range) do concentration
		hill(concentration,tf.K,tf.n)
	end
end

#-----------------------------------------------------------------------------------------
"""
	animate_hill()

Create a graphical animation of the Hill function as the value of the cooperation (n)
changes from 1 up to 10.
"""
function animate_hill()
	tf_range = 0:30
	K = 10

	anim_n = Observable(1.0)           		# The animated cooperation value is Observable.
	hill_curve = map(anim_n) do n			# Create an expression curve that listens
		expression(HillTF(tf_range,K,n))    # directly to the current value (n) of anim_n.
	end

	fig = lines(tf_range, hill_curve)       # Create a plot of the (current) hill_curve ...
	display(fig)                        	# ... and display the result

	for n in 1:0.1:10                   	# Finally, animate the figure: Step through the
		anim_n[] = n                   		# values of n, and watch the hill_curve change
		sleep(0.1)                      	# as the value of anim_n moves from 1 to 10.
	end
end

#-----------------------------------------------------------------------------------------
"""
	fun_graphics( permutation::AbstractVector{Int}=ones(5), evaluate=false)

Print out the given permutation of a particular block of scrambled graphics code, then evaluate
that permuted code if requested.
"""
function fun_graphics( permute::AbstractVector{Int}=1:11; evaluate=false)
	if length(permute) != 11 || any(map(x->!(x in 1:11), permute))
		# Bad permutation:
		return false
	end

	permuted_code = scrambled_code[permute]
	for line in permuted_code
		println(line)
	end
	println()

	if [sum(permute[3:5]),permute[[1,2,6,7,11]]...,sum(permute[8:10])] != [17,9,2,3,8,7,20]
		# Permuted code is not executable:
		return false
	end

	if evaluate
		display( eval(Meta.parse(join(permuted_code))))
	end

	true
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Demonstrate the use of the HillTF module.
"""
function demo()
	tf = HillTF(0:30,5)
	expression(tf)
end

end		# of FunGraphics