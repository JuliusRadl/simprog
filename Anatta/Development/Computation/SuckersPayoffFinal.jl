#==============================================================================#
"""
	SuckersPayoff

Module SuckersPayoff: A simulation of different strategies interacting under
varied payoff scenarios.

This program investigates the following research question:
	Does the price for cooperation in the face of defection affect how reliably
	cooperative strategies can prevail?

Authors: Sophia Prögler, Julius Radl, 2026-01-18
"""
module SuckersPayoff

using Random, CairoMakie

#-------------------------------------------------------------------------------
# Module types:
#-------------------------------------------------------------------------------
"""
	SimulationModel

Contains all necessary parameters for a simulation run via an iterated prisoners
dilemma. Holds simulation results.
"""
mutable struct SimulationModel
	payoffs::Matrix{Float64}		# Payoff matrix
	strategies::Matrix{Float64}		# Strategies of current population
	x::Vector{Float64}				# Current frequency of strategies
	novak::Matrix{Float64}			# Novak's matrix of long term payoffs
	nTrials::Int64					# Number of Trials over which we will average results
	nMutations::Int64				# How often bad strategies will be replaced
	nGenerations::Int64				# Simulation steps between mutation events
	mu::Float64						# Fraction of strategies to be replaced on mutation event
	initStrat::Vector{Float64}		# Initial strategy seed
	initEpsilon::Float64			# Initial strategy variation around seed
	history::Matrix{Float64}			# Simulation results averaged over trials

	"""
		SimulationModel(
			payoffs,
			nStrategies,
			nTrials,
			nMutations,
			nGenerations,
			mu,
			initStrat,
			initEpsilon
		)

	Constructor for SimulationModel object with sensible default values.
	Reserves space and evenly distributes population frequency.
	"""
 	function SimulationModel(
		payoffs = [4 0; 5 1],
		nStrategies = 20,
		nTrials = 20,
		nMutations = 100,
		nGenerations = 20,
		mu = 0.1,
		initStrat = [0.075, 0.075],
		initEpsilon = 0.025)

		new(
			Float64.(payoffs),
			zeros(Float64, nStrategies, 2),			# Empty p,q-Matrix for nStrategies
			ones(nStrategies)/nStrategies,			# Uniform strategy distribution
			zeros(Float64,2,2),						# novak matrix placeholder
			nTrials,
			nMutations,
			nGenerations,
			mu,
			initStrat,
			initEpsilon,
			zeros(Float64, nMutations+1, 2)			# Empty history matrix
		)
	end
end

#-------------------------------------------------------------------------------
# Module methods:
#-------------------------------------------------------------------------------
"""
	mutate!(sm)

Replace the least fit fraction mu of the current strategies by new, randomly
selected strategies and adjust their frequency to give them a chance.
"""
function mutate!(sm::SimulationModel)

	# Number of worst Strategies
	nWorstStrats = ceil(Integer, length(sm.x) * sm.mu)

	# Find worst and best strategies
	ranking = sortperm(sm.x)
	worstStratsIndex = ranking[1:nWorstStrats]
	bestStratsIndex = ranking[nWorstStrats+1:end]

	# Generate and introduce new strategies
	sm.strategies[worstStratsIndex, :] = rand(Float64, nWorstStrats, 2)

    # Give new strategies an even frequency share, adjust overall frequency
	sm.x[worstStratsIndex] .= (1/length(sm.x))
	remaining_freq = 1 - sum(sm.x[worstStratsIndex])
	sm.x[bestStratsIndex] .= (
		remaining_freq *
		(sm.x[bestStratsIndex] / sum(sm.x[bestStratsIndex]))
	)

	# Recalculate Payoff matrix
	novak!(sm)
end

#-------------------------------------------------------------------------------
"""
	setSuckers!(sm::simulationModel)

Helper function for adjusting the sucker's payoff. Increases S by a fraction of
the interval between S and P (Punishment for mutual defection). The Fraction
must not be set to 1 to maintain PD prerequisites.
"""
function setSuckers!(sm::SimulationModel, fraction)

	sm.payoffs[1, 2] += (sm.payoffs[2, 2] - sm.payoffs[1, 2]) * fraction
end

#-------------------------------------------------------------------------------
"""
	set!(sm::SimulationModel)

Set all current strategies to the given strategy with random variation epsilon,
then evenly distribute population frequency among strategies, making sure that
p,q are within the interval ]0,1[.
"""
function set!(sm::SimulationModel)
	
	nStrategies = length(sm.x)
	strategies = (
		(2 * rand(Float64, nStrategies, 2) .- 1) *
		sm.initEpsilon .+ sm.initStrat'
	)
	sm.strategies = clamp.(strategies, eps(), 1-eps())
	sm.x .= 1/nStrategies
end

#-------------------------------------------------------------------------------
"""
	simulate!(sm::SimulationModel)

Simulate the effect of repeated interaction for nGenerations steps using
Runge-Kutta-2 and the replicator equation. Update strategy frequencies
accordingly.
"""
function simulate!(sm::SimulationModel)
	
	dt = 1.0
	for _ in 1:sm.nGenerations
		# Perform RK2 half-step:
		R = sm.x' * sm.novak * sm.x
		dxdt = sm.x .* (sm.novak * sm.x .- R)
		half_step = sm.x + dxdt * (dt/2)
		
		# Perform RK2 full-step:
		new_R = half_step' * sm.novak * half_step
		new_dxdt = half_step .* (sm.novak * half_step .- new_R)
		sm.x = sm.x + new_dxdt * dt

		# Normalize frequencies
		sm.x = sm.x ./ sum(sm.x)
	end
end

#-------------------------------------------------------------------------------
"""
	novak!(sm::SimulationModel)

Recalculate Novak's expected payoff from current strategies and payoff matrix.
"""
function novak!(sm::SimulationModel)

	# Legibility
	Te = sm.payoffs[2, 1]						# Temptation to defect
	Re = sm.payoffs[1, 1]						# Reward for mutual cooperation
	Pu = sm.payoffs[2, 2]						# Punishment for mutual defection
	Su = sm.payoffs[1, 2]						# Sucker's payoff
	p = sm.strategies[:, 1]						# Reciprocation
	q = sm.strategies[:, 2]						# Forgiveness
	
	r = p .- q
	s_ij = (r .* q' .+ q) ./ (1 .- (r * r'))	
	s_ji = s_ij'								
	
	Te_part = Te * ((1 .- s_ij) .* s_ji)
	Re_part = Re * (s_ij .* s_ji)
	Pu_part = Pu * ((1 .- s_ij) .* (1 .- s_ji))
	Su_part = Su * (s_ij .* (1 .- s_ji))

	sm.novak = Te_part + Re_part + Pu_part + Su_part
end

#-------------------------------------------------------------------------------
"""
	avgstrategy(sm::SimulationModel)

Return vector containing average strategy across current population.
"""
function avgstrategy(sm::SimulationModel)
	(sm.x'*sm.strategies)'
end

#-------------------------------------------------------------------------------
"""
	resetHistory!(sm::SimulationModel)

Set all values in history to 0.0.
"""
function resetHistory!(sm::SimulationModel)
	sm.history = zeros(Float64, sm.nMutations+1, 2)
end

#-------------------------------------------------------------------------------
"""
	run!(sm::SimulationModel)

Explore which strategies survive. Run nTrials and average results over trials.
"""
function run!(sm::SimulationModel)

	resetHistory!(sm)

	for i in 1:sm.nTrials
		print("[Trial $(i)/$(sm.nTrials)]\r")
		set!(sm)									# Reset to initial exploitative population
		sm.history[1, :] += avgstrategy(sm)			# Save initial strategy average to history
		
		for mut in 1:sm.nMutations
			mutate!(sm)	
			simulate!(sm)							# Iterate PD interactions
			sm.history[mut+1, :] += avgstrategy(sm)	# Accumulate new simulation results
		end
	end
	sm.history = sm.history ./ sm.nTrials			# Normalize results
end

#-------------------------------------------------------------------------------
# Display methods:
#-------------------------------------------------------------------------------
"""
	present(models::Vector)

Present the results of our investigation.
"""
function present(models::Vector{SimulationModel})

	fig = Figure(
		fontsize=18,
		regular="Helvetica",
		linewidth=5,
		size=(1412, 1000)
	)

	present_title!(fig)
	present_context!(fig)
	present_problem!(fig)
	present_method!(fig)
	present_results!(fig)
	present_plots!(fig, models, 4)
	present_implications!(fig)
	present_references!(fig)
	resize_to_layout!(fig)

	path = joinpath(pwd(), "SuckersPayoffInfoSheet.pdf")
	save(path, fig)
	println("[!] Saved poster to $(path)")
end

#-------------------------------------------------------------------------------
"""
	present_title!(fig)

In the given figure, show title, subtitle and authors.
"""
function present_title!(fig)

	Label(
		fig[1, 1][1, 1:4],
		fontsize = 30,
		justification=:left,
		halign=:left,
		color=:black,
		font=:bold,
		tellwidth = false,
		text = "Hängt die Entwicklung von Kooperation vom\nSucker's Payoff ab?"
	)

	Label(
		fig[1, 1][2, 1:4],
		fontsize = 20,
		justification=:left,
		halign=:left,
		color=:black,
		font=:italic,
		tellwidth = false,
		text = "Computational Game Dynamics angewendet auf " *
			"Evolutionsbiologie.\nSophia Prögler & Julius Radl, Hochschule " *
			"Weihenstephan-Triesdorf, 18.01.2026"
	)
end

#--------------------------------------------------------------------------------
"""
	present_context!(fig)

In the given figure, describe the relevance of the Cooperation model in the
context of current scientific debate.
"""
function present_context!(fig)

	Box(
		fig[1, 1][1:2, 5:7],
		cornerradius = (0, 0, 10, 0),
		color = Makie.cgrad(:cool, alpha = 0.2)[12],
		strokecolor=:black,
		strokewidth=2
	)

	Label(
		fig[1, 1][1:2, 5:7],
		justification=:left,
		halign=:left,
		padding = 10,
		color=:black,
		word_wrap = true,
		text = "Kontext:\n" *
			"Kooperative Strategien können sich unter geeigneten Bedingungen " *
			"in einer von Verrätern dominierten Population langfristig " *
			"durchsetzen (Nowak, 2006). Den Erfolg solcher Strategien " *
			"können wir in der Natur beobachten."
	)
end

#-------------------------------------------------------------------------------
"""
	present_problem!(fig)

In the given figure, describe the specific problem that the Cooperation model
addresses in current scientific debate.
"""
function present_problem!(fig)

	Box(
		fig[2, 1][1, 4:5],
		cornerradius = (0, 0, 10, 10),
		color = Makie.cgrad(:cool, alpha = 0.2)[30],
		strokecolor=:black,
		strokewidth=2,
	)
	Label(
		fig[2, 1][1, 4:5],
		justification=:left,
		halign=:left,
		padding=10,
		color=:black,
		word_wrap = true,
		text = "Problem:\n" *
			"Welchen Einfluss nimmt der Preis, den Kooperateure bei der " *
			"Interaktion mit Verrätern zahlen, auf ihre Fähigkeit, sich " *
			"langfristig durchzusetzen? Diesen Preis können wir durch den " *
			"Sucker’s Payoff (S) aus der evolutionären Spieltheorie " *
			"abbilden: Ein höheres S bedeutet einen niedrigeren Preis. " *
			"Dabei nehmen wir an, dass sich kooperative Strategien schneller " *
			"durchsetzen können, je höher S ist."
	)
end

#-------------------------------------------------------------------------------
"""
	present_method!(fig)

In the given figure, describe the exact, algorithmic sequence of steps and
measurements that others scientists must use in order to reproduce our findings.
"""
function present_method!(fig)

	Box(
		fig[2, 1][1, 1:3],
		cornerradius = (10, 10, 0, 0),
		color = Makie.cgrad(:cool, alpha = 0.2)[50],
		strokecolor=:black,
		strokewidth=2
	)
	Label(
		fig[2, 1][1, 1:3],
		justification=:left,
		halign=:left,
		padding=10,
		color=:black,
		word_wrap = true,
		text = "Methode:\n" *
			"Typenfitness wird in Populationen durch Interaktion bestimmt " *
			"(Hofbauer & Sigmund, 1998). Wir simulieren rechnergestützt " *
			"(Nowak 2006) eine mutierende Population im iterierten " *
			"Gefangenendilemma (PD), in der Interaktionsstrategien durch " *
			"Wahrscheinlichkeiten für Reziprozität (p) und Vergebung (q) " *
			"definiert sind. In der Auszahlungsmatrix des PD gilt: " *
			"Einseitige Kooperation (S) darf sich nicht mehr lohnen als " *
			"beidseitiger Verrat (P) (Axelrod 1984). Deshalb nähern wir S " *
			"schrittweise an P an. Wir betrachten die über viele Episoden " *
			"gemittelte Entwicklung der durchschnittlichen p- & q-Werte und " *
			" diskutieren anschließend den Einfluss von S auf die Frequenz " *
			"kooperativer Strategien."
	)
end


#-------------------------------------------------------------------------------
"""
	present_results!(fig)

In the given figure, describe the precise, uninterpreted and error-prone
empirical data that arose from implementing the method.
"""
function present_results!(fig)

	Box(
		fig[3, 1][1, 1][1, 1],
		cornerradius = (10, 10, 0, 0),
		color = Makie.cgrad(:cool, alpha = 0.2)[75],
		strokecolor=:black,
		strokewidth=2
	)
	Label(
		fig[3, 1][1, 1][1, 1],
		justification=:left,
		halign=:left,
		padding=10,
		color=:black,
		word_wrap = true,
		text = "Ergebnisse:\n" *
			"Solange S unter etwa der Hälfte von P liegt, fällt q anfangs " *
			"ab, während p steigt, dann holt q aber auf. Sobald S darüber " *
			"liegt, steigt q von Anfang an. In allen Graphen lässt sich der " *
			"gleiche Trend erkennen: Früher oder später wird p = 1 und " *
			"q = 0.5 erreicht. Je höher S ist, auf desto weniger Umwegen " *
			"erreicht die Population diesen Zustand."
	)
end

#-------------------------------------------------------------------------------
"""
	present_implications!(fig)

In the given figure, describe the meaning and implications of our results in
relation to the previous description of the problem and its context.
"""
function present_implications!( fig)
	Box(
		fig[3, 1][1, 1][2, 1],
		cornerradius = (10, 10, 0, 0),
		color = Makie.cgrad(:cool, alpha = 0.2)[100],
		strokecolor=:black,
		strokewidth=2
	)
	Label(
		fig[3, 1][1, 1][2, 1],
		justification=:left,
		halign=:left,
		padding=10,
		color=:black,
		word_wrap = true,
		text = "Implikationen:\n" *
			"Wir schließen daraus, dass der Preis für einseitige Kooperation " *
			"einen entscheidenden Einfluss auf die Evolution nimmt: Je näher " *
			"er am Preis für beidseitigen Verrat liegt, desto zuverlässiger " *
			"und schneller setzen sich kooperative Strategien durch. In " *
			"einer solchen Umgebung wird der Hang zu Vergebung weniger hart " *
			"bestraft und Reziprozität genährt."
	)
end

#-------------------------------------------------------------------------------
"""
	plot_history!(fig, sm::SimulationModel)

Create axis for the average strategy history of a single SimulationModel.
"""
function plot_history!(fig, sm::SimulationModel)

	axis = Axis(
		fig,
		xgridvisible = true,
		ygridvisible = true,
		xticks = 0:0.2:1,
		yticks = 0:0.2:1,
		xticklabelsize = 12,
		yticklabelsize = 12,
		aspect = DataAspect()
	)

	lines!(
		axis,
		sm.history[:, 1],
		sm.history[:, 2],
		color = LinRange(1, 10, length(sm.history[:, 1])),
		colormap=Reverse(:copper),
		linewidth = 5
	)

	xlims!(axis, 0, 1)
	ylims!(axis, 0, 1)

	text!(
		axis,
		0.05, 0.95,
		text = "S = $(round(sm.payoffs[1,2], digits = 2))",
		align = (:left, :top),
		space = :relative,
		fontsize = 14,
		color = :black,
	)
end

#-------------------------------------------------------------------------------
"""
	present_plots!(fig, models::Vector{SimulationModel}, nCols = 4)

Plot histories of multiple simulation models and arrange them in a grid with
nCols columns, a title, a colorbar and a description of dimensions.
"""
function present_plots!(fig, models::Vector{SimulationModel}, nCols = 4)

	Label(
		fig[3, 1][1, 2:4][1, 1],
		color=:black,
		font=:bold,
		tellwidth=false,
		text = "Entwicklung von Reziprozität und Vergebung bei P = 1.0"
	)

	Colorbar(
		fig[3, 1][1, 2:4][2, 1][1:2, 1],
		limits = (0, models[1].nMutations),
		label = "Zeit",
		colormap = Reverse(:copper),
		flipaxis = false,
		ticklabelsvisible = false
	)
	
	# Distribute plots across columns
	nModels = length(models)
	for i in 1:nModels
		row = div(i-1, nCols) + 1
		col = mod(i-1, nCols) + 2
		plot_history!(fig[3, 1][1, 2:4][2, 1][row, col], models[i])
	end

	Label(
		fig[3, 1][1, 2:4][3, 1],
		color=:black,
		tellwidth = false,
		padding = (0, 0, 0, -15), 
		text = "[x] = Reziprozität, [y] = Vergebung"
	)
end

#-------------------------------------------------------------------------------
"""
	present_references!(fig)

In the given figure, list all sources of information cited in all the above
descriptions.
"""
function present_references!( fig)

	Box(
		fig[4, 1],
		color = Makie.cgrad(:cool, alpha = 0.03)[100],
		strokecolor=:black,
		strokewidth=2,
		tellwidth = false
	)

	Label(
		fig[4, 1],
		justification=:left,
		halign=:left,
		padding=10,
		color=:black,
		word_wrap = true,
		fontsize = 11,
		tellwidth = false,
		text="Quellen: " *
			"Axelrod, R. (1984) The evolution of cooperation. " *
			"Basic Books / " *
			"Hofbauer, J. & K. Sigmund (1998). Evolutionary Games and " *
			"Population Dynamics. Cambridge University Press / " *
			"Nowak, M. A. (2006) Evolutionary dynamics. " *
			"Harvard University Press"
	)
end

#-------------------------------------------------------------------------------
# Demo methods:
#-------------------------------------------------------------------------------
"""
	demo()

Use PD dynamics to investigate the following research question:
	Does the price for cooperation in the face of defection affect how reliably
	cooperative strategies can prevail?
"""
function demo()

	# Arguments for SimulationModel constructor
	payoffs = [4 0; 5 1]
	nStrategies = 50
	nTrials = 20
	nMutations = 300
	nGenerations = 50
	mu = 0.01		
	initStrat = [0.075, 0.075]
	initEpsilon = 0.025

	# Suckers' Payoff variation
	sChange = collect(range(0.0, 0.99, 8))

	nRuns = length(sChange)
	models = Vector{SimulationModel}(undef, nRuns)
	for run in 1:nRuns
		models[run] = SimulationModel(
			payoffs,
			nStrategies,
			nTrials,
			nMutations,
			nGenerations,
			mu,
			initStrat,
			initEpsilon
		)
		setSuckers!(models[run], sChange[run])
		println(
			"[+] Now running Simulation for model $(run)/$(nRuns): " *
			"S = $(round(models[run].payoffs[1, 2], digits = 2))"
		)			
		run!(models[run])
		println("[+] Finished Simulation for model $(run)/$(nRuns)")
	end
	println("[!] Finished simulating all models")
	present(models)
end

end		# ... of module SuckersPayoff