#========================================================================================#
"""
	Selectors

Module Selectors: A model of super- and sublinear selection in a replicating population.

Author: Julius Radl

Palfreymans Lösung:

Hat plot3!(axis, selector)-Methode IN Selector, in der
Simplex.plot3!(axis) und Simplex.plot3!(axis, selector)
aufgerufen wird. Simplex.plot3!(axis) zeichnet NUR das
Simplex-Dreieck, nichts weiter. Dh ich muss immer beide
aufgerufen.

Meine Fehler bei Runge-Kutta-2: 
-hab dxdt nicht halbiert für die Berechnung des Halbschritts
-hab R nur einmal berechnet, aber beim ganzen Schritt muss ich
R neu aus dem Ergebnis des Halbschritts berechnen!
"""
module Selectors

include( "Simplex.jl")

using GLMakie, LinearAlgebra
#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Selector

A Selector represents the time-evolution of a sub- or super-linear selection model.

Palfreymans Lösung: 

Hat noch 3 Variablen:
-epsilon (wird von 1 abgezogen / addiert => damit c beschrieben)
-resolution (Zahl der Zeitschritte im Zeitbereich)
-timescale (Dauer als eigener Vektor)

Initialisiert alle Variablen (auch freqs). Dadurch muss der Vektor nicht nachher
bei jedem Zeitschritt erweitert werden, sondern wird nur gefüllt (braucht also
kein push!)

Normiert spezifische Wachstumsraten

Hat einen default-Wert für epsilon

timescale ist ein Vektor, der als resolution+1 initialisiert wird
(Brauche ich nicht für den Simplex, aber wenn ich später einen Graphen
mit Zeitachse erstellen wollen würde, dann wären passende Werte gut)
"""
struct Selector
	sgr::Vector{Float64}					# Specific growth rates of three types in min^-1
	c::Float64								# Exponent determining type of linearity
	freqs::Vector{Vector{Float64}}			# Evolution of the frequencies of each type

	function Selector(sgr::Vector{Float64}, c)
		freqs = [[0.0, 0.0, 0.0]]			# Initialize variable
		new(sgr, c, freqs)
	end
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------

"""
	simulate!(sel, init_freq, time_steps)

Simulates evolution of populations in sel with initial frequencies init_freq
over the duration time_steps using Runge-Kutta-2 and saves each step to sel.
"""

function simulate!(sel::Selector, init_freq::Vector{Float64}, time_steps::Int64)
	empty!(sel.freqs)

	push!(sel.freqs, init_freq)

	dt = 1.0

	for _ in 1:time_steps

		# für Lesbarkeit
		current_freqs = sel.freqs[end]

		# R (durchschnittliche Wachstumsrate) berechnen
		R = dot(sel.sgr, current_freqs.^sel.c)
		
		# Wachstumsraten berechnen pro Typ & Runge Kutta 2
		dxdt = (sel.sgr .* current_freqs) .- (R * current_freqs)
		new_freqs = rktwo_step(dt, dxdt, current_freqs)

		# Aktuelle Frequenzen normalisieren
		new_freqs = new_freqs ./ sum(new_freqs)

		# Aktuelle Frequenzen eintragen
		push!(sel.freqs, new_freqs)
	end
end

#-----------------------------------------------------------------------------------------
"""
	rktwo_step(dt::Float64, dxdt::Float64, x0::Vector{Float64})

Calculates a single step of the Runge-Kutta-2-Method for values
in a vector.
"""
function rktwo_step(dt::Float64, dxdt::Vector{Float64}, x0::Vector{Float64})


		half_step = x0 .+ x0 .* dxdt * (dt/2)
		x1 = x0 .+ half_step .* dxdt * dt

		return x1
end

#-----------------------------------------------------------------------------------------
"""
	unittest()

Unit-test the Selectors module.
"""
function unittest()
	println("\n============ Unit test Selectors: ===============")

	println("First we will need some Figures and Axes to show our plots off:")
	fig = Figure()
    ax1 = Axis(fig[1, 1])
	ax2 = Axis(fig[1, 2])
	ax3 = Axis(fig[2, 1])
	ax4 = Axis(fig[2, 2])
	display(fig)
	
	println("Let's start with a linear selection model (c = 1).",
	"First, we will need a Selector, defining its growth rates and c-factor:")
	linear = Selector([0.4, 0.5, 0.6], 1.0)
	display(linear)

	println("Now we will simulate and plot the evolution of its types' populations,",
	"using the their initial frequencies and a number of timesteps:")
	simulate!(linear, [0.6, 0.3, 0.1], 1000)

	println("Let's show our results off:")
	Simplex.plot3!(ax1)
	Simplex.plot3!(ax1, linear.freqs)

	println("A sublinear Model (c=0.5): Only the specific growth rates",
	"determine where the equilibrium point ends up!")
	sublinear = Selector([0.2, 0.4, 0.5], 0.5)
	display(sublinear)
	simulate!(sublinear, [0.6, 0.3, 0.1], 1000)
	Simplex.plot3!(ax2)
	Simplex.plot3!(ax2, sublinear.freqs)

	println("A sublinear Model (c=0.5), but with different initial populations")
	sublinear2 = Selector([0.2, 0.4, 0.5], 0.5)
	display(sublinear2)
	simulate!(sublinear2, [0.1, 0.3, 0.6], 1000)
	Simplex.plot3!(ax3)
	Simplex.plot3!(ax3, sublinear2.freqs)

	println("A superlinear Model (c=1.5): The biggest initial population wins!")
	superlinear = Selector([0.4, 0.5, 0.6], 1.5)
	display(superlinear)
	simulate!(superlinear, [0.6, 0.3, 0.1], 1000)
	Simplex.plot3!(ax4)
	Simplex.plot3!(ax4, superlinear.freqs).parent
end

end		# ... of module Selectors