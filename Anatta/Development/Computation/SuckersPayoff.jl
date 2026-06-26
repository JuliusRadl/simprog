#========================================================================================#
"""
	SuckersPayoff

Module SuckersPayoff: A simulation of different strategies interacting under under varied
payoff scenarios.

This program investigates the following research question:
	Can cooperative strategies still prevail, even when the price for cooperation in the face
	of defection increases?

Authors: Sophia Prögler, Julius Radl, 2026-01-18
"""
module SuckersPayoff

using Random, CairoMakie, LinearAlgebra

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	SimulationModel

Contains all necessary parameters for a simulation run via an iterated prisoners dilemma.
Holds simulation results.
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
	history::Vector{Vector}			# Simulation results averaged over trials

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
			Float64.(payoffs),						# Payoff matrix
			zeros(Float64, nStrategies, 2),			# Empty p,q-Matrix for nStrategies
			ones(nStrategies)/nStrategies,			# Uniform strategy distribution
			zeros(Float64,2,2),						# novak matrix placeholder
			nTrials,
			nMutations,
			nGenerations,
			mu,
			initStrat,
			initEpsilon,
			Vector{Vector}(undef, nMutations+1)		# History placeholder
		)
	end
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	mutate!(sm)

Replace the least fit fraction mu of the current strategies by new, randomly selected
strategies and adjust their frequency to give them a chance.
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

    # Give new strategies an even frequency share and andjust
	# old stratgies' frequencies
	sm.x[worstStratsIndex] .= (1/length(sm.x))
	remaining_freq = 1 - sum(sm.x[worstStratsIndex])
	sm.x[bestStratsIndex] .= remaining_freq * (
		sm.x[bestStratsIndex] / sum(sm.x[bestStratsIndex]))

	# Recalculate Payoff matrix
	novak!(sm)
end

#-----------------------------------------------------------------------------------------
"""
	setS!(sm::simulationModel)

Helper function for adjusting the sucker's payof.
"""
function setS!(sm::SimulationModel, S)

	sm.payoffs[1, 2] = S
end

#-----------------------------------------------------------------------------------------
"""
	set!(sm::SimulationModel)

Set all current strategies to the given strategy with random variation epsilon,
then evenly distribute population frequency among strategies, making sure that
p,q are within the interval ]0,1[.
"""
function set!(sm::SimulationModel)
	
	nStrategies = length(sm.x)
	strategies = (2 * rand(Float64, nStrategies, 2) .- 1) * sm.initEpsilon .+ sm.initStrat'
	sm.strategies = clamp.(strategies, eps(), 1-eps())
	sm.x .= 1/nStrategies
end

#-----------------------------------------------------------------------------------------
"""
	simulate!(sm::SimulationModel)

Simulate the effect of repeated interaction for nGeneration steps.
"""
function simulate!(sm::SimulationModel)

	# for _ in 1:sm.nGenerations
	# 	fitness = sm.novak * sm.x
	# 	average_fitness = sm.x' * fitness
	# 	sm.x .= sm.x .*(fitness ./ average_fitness)
	# end

	dt = 1.0
	for _ in 1:sm.nGenerations

		# Perform RK2 half-step:
		R = sm.x' * sm.novak * sm.x
		dxdt = sm.x .* (sm.novak * sm.x .- R)
		half_step = sm.x + sm.x .* dxdt * (dt/2)
		
		# Perform RK2 full-step:
		new_R = half_step' * sm.novak * half_step
		new_dxdt = half_step .* (sm.novak * half_step .- new_R)
		x = sm.x + sm.x .* new_dxdt * dt

		# Normalize frequencies
		sm.x = sm.x ./ sum(sm.x)
	end

end

#-----------------------------------------------------------------------------------------
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

	r = p - q
	s_ij = (r .* q' .+ q) ./ (1 .- r * r')								

	Re_part = Re * (s_ij .* s_ij')
	Su_part = Su * (s_ij .* (1 .- s_ij'))
	Te_part = Te * ((1 .- s_ij) .* s_ij')
	Pu_part = Pu * ((1 .- s_ij) .* (1 .- s_ij'))

	sm.novak = Te_part + Re_part + Pu_part + Su_part
end

#-----------------------------------------------------------------------------------------
"""
	avgstrategy(sm::SimulationModel)

Return vector containing average strategy across current population.
"""
function avgstrategy(sm::SimulationModel)
	(sm.x'*sm.strategies)'
end

#-----------------------------------------------------------------------------------------
"""
	resetHistory!(sm::SimulationModel)

Set all values in history to 0.0.
"""
function resetHistory!(sm::SimulationModel)
	sm.history = [[0.0,0.0] for _ in 1:sm.nMutations+1]
end

#-----------------------------------------------------------------------------------------
"""
	run!(sm::SimulationModel)

Explore which strategies survive. Run nTrials and average results over trials.
"""
function run!(sm::SimulationModel)

	resetHistory!(sm)

	for _ in 1:sm.nTrials
		set!(sm)									# Reset to initial exploitative population
		sm.history[1] += avgstrategy(sm)			# Save initial strategy average to history

		for mut in 1:sm.nMutations
			mutate!(sm)								# Replace current weakest fraction of strategies
			simulate!(sm)							# Iterate PD interactions
			sm.history[mut+1] += avgstrategy(sm)	# Accumulate new simulation results
		end
	end
	sm.history = sm.history ./ sm.nTrials			# Normalize results
end

#-----------------------------------------------------------------------------------------
# Display methods:
#-----------------------------------------------------------------------------------------
# TODO brauchen wir das wirklich?
"""
	behaviour(p::Float64, q::Float64)

Return 4-char string classifying the strategy S(p,q) according to cooperating/defecting
behavioural type:

	(0.0,0.0) : AllD - Always defect
	(0.0,0.5) : Xplt - Exploiting
	(0.0,1.0) : Ctry - Contrary
	(0.5,0.0) : STFT - Selfish TFT
	(0.5,0.5) : Gmbl - Gambling
	(0.5,1.0) : Plac - Placating
	(1.0,0.0) : RTFT - Rulebook Tit-for-tat
	(1.0,0.5) : GTFT - Generous TFT
	(1.0,1.0) : AllC - Always cooperate
"""
function behaviour(strategy::Vector)
	p,q = strategy
	cutoff = 0.2
	
	if q > 1-cutoff
		return (p<cutoff) ? "Ctry" : ((p<1-cutoff) ? "Plac" : "AllC")
	elseif q > cutoff
		return (p<cutoff) ? "Xplt" : ((p<1-cutoff) ? "Gmbl" : "GTFT")
	else
		return (p<cutoff) ? "AllD" : ((p<1-cutoff) ? "STFT" : "RTFT")
	end
end

#-----------------------------------------------------------------------------------------
"""
	show(sm::SimulationModel, bpause)

Display the current Players frequencies at the console in descending frequency order.
"""
function show(sm::SimulationModel, proportion=1, bpause::Bool=false)
	proportion = max(0,min(1,proportion))
	ntoshow = ceil(Int,proportion*length(sm.x))
	toshow = sortperm(sm.x,rev=true)[1:ntoshow]
	for i in 1:ntoshow
		println("$(behaviour(sm.strategies[toshow[i], :])): $(sm.x[toshow[i]])")
	end
	println()
	if bpause
		# Pause to allow user to digest output
		readline()
	end
end

#-----------------------------------------------------------------------------------------
"""
	present(models::Vector)

Present the results of our investigation.
"""
function present(models::Vector)
	fig = Figure(fontsize=18, regular="Helvetica", linewidth=5, size=(1412, 1000))

	present_title!(fig)
	present_context!(fig)
	present_problem!(fig)
	present_method!(fig)
	present_results!(fig)
	present_plots!(fig, models, 4)
	present_implications!(fig)
	present_references!(fig)

	resize_to_layout!(fig)

	fig
	# TODO: später reinkommentieren
	# save("CooperationInfoSheet.pdf",fig)
	# nothing
end

#-----------------------------------------------------------------------------------------
"""
	present_title!(fig)

In the given figure, show title, subtitle and authors.
"""
function present_title!(fig)
	Label(fig[1, 1][1, 1:5], fontsize = 30, justification=:left, halign=:left, color=:black, font=:bold,
		text = "Hängt die Entwicklung von Kooperation vom Sucker's Payoff ab?"
	)
	Label(fig[1, 1][2, 1:5], fontsize = 20, justification=:left, halign=:left, color=:black, font=:italic,
		text = "Computational Game Dynamics angewendet auf Evolutionsbiologie.\n" *
			"Sophia Prögler, Julius Radl, Hochschule Weihenstephan-Triesdorf, 13.01.2026"
	)
end

#-----------------------------------------------------------------------------------------
"""
	present_context!(fig)

In the given figure, describe the relevance of the Cooperation model in the context of current
scientific debate.
"""
function present_context!(fig)
	Box(fig[1, 1][1:2, 6:7], cornerradius = (0, 0, 10, 0), color = Makie.cgrad(:cool, alpha = 0.2)[10],
		strokecolor=:black, strokewidth=2)
	Label(fig[1, 1][1:2, 6:7], justification=:left, padding = 10, halign=:left, color=:black,
		word_wrap = true,
		text = "Context:\n" *
			"Bereits 1871 formulierte Charles Darwin in The Descent of Man " *
			"ein fundamentales Paradoxon des sozialen Zusammenlebens."
	)
end

#-----------------------------------------------------------------------------------------
"""
	present_problem!( fig)

In the given figure, describe the specific problem that the Cooperation model addresses in current
scientific debate.
"""
function present_problem!(fig)

	Box(fig[2, 1][1, 1], cornerradius = (10, 10, 0, 0), color = Makie.cgrad(:cool, alpha = 0.2)[25],
		strokecolor=:black, strokewidth=2, tellheight = false)
	Label(fig[2, 1][1, 1], justification=:left, halign=:left, padding=10, color=:black, word_wrap = true,
		tellheight = false,
		text = "Problem:\n" *
			"Kooperative Strategien können sich unter gegebenen Payoffs durchsetzen. " *
			"Wie ist diese Entwicklung von der Ausprägung des S-Werts abhängig?"

	)
end

#-----------------------------------------------------------------------------------------
"""
	present_method!(fig)

In the given figure, describe the exact, algorithmic sequence of steps and measurements that others
scientists must use in order to reproduce our findings.
"""
function present_method!(fig)
	Box(fig[2, 1][1, 2:4], cornerradius = (0, 0, 10, 10), color = Makie.cgrad(:cool, alpha = 0.2)[40],
		strokecolor=:black, strokewidth=2)
	Label(fig[2, 1][1, 2:4], justification=:left, halign=:left, padding=10, color=:black, word_wrap = true,
		text = "Method:\n" *
			"Die Untersuchung setzt voraus, dass Organismen aktive, konstruktive Rollen" *
			"in ihrer eigenen Entwicklung einnehmen (Laland 2024), wodurch extragenetische " *
			"Prozesse die Richtung der Evolution beeinflussen können. Wir simulieren " *
			"rechnergestützt (Nowak 2006) eine mutierende Population im iterierten " *
			"Gefangenendilemma (PD), in der Strategien durch ihre Wahrscheinlichkeiten " *
			"für Kooperation (p) und Vergebung (q) definiert sind. Um die Stabilität " *
			"dieses Systems zu prüfen, verändern wir graduell die Auszahlungsmatrix: " *
			"Wir variieren systematisch den Suckers' Payoff (S) und beobachten den Effekt " *
			"auf den Erfolg kooperativer Strategien."
	)
end


#-----------------------------------------------------------------------------------------
"""
	present_results!( fig)

In the given figure, describe the precise, uninterpreted and error-prone empirical data that
arose from implementing the method.
"""
function present_results!(fig)
	Box(fig[3, 1][1, 1][1, 1], cornerradius = (10, 10, 0, 0), color = Makie.cgrad(:cool, alpha = 0.2)[55],
		strokecolor=:black, strokewidth=2)
	Label(fig[3, 1][1, 1][1, 1], justification=:left, halign=:left, padding=10, color=:black, word_wrap = true,
	text = "Ergebnisse:\n" *
		"Bei S >= 0 steigen die durchschnittlichen p- und q-Werte gleichmäßig " *
		"an, bis q ein Maximum von etwa 0.5 erreicht. Im Fall S = 0 wird nur " *
		"ein Maximum von q ≈ 0.35 erreicht. Ab da fallen die q-Werte wieder, " *
		"während p weiter zu steigen scheint. Bei negativem S steigen p und q nur " *
		"noch geringfügig, bevor beide wieder fallen. Der Effekt ist bei geringerem S ausgeprägter."
	)
end

#-----------------------------------------------------------------------------------------
"""
	present_implications!( fig)

In the given figure, describe the meaning and implications of our results in relation to the
previous description of the problem and its context.
"""
function present_implications!( fig)
	Box(fig[3, 1][1, 1][2, 1], cornerradius = (10, 10, 0, 0), color = Makie.cgrad(:cool, alpha = 0.2)[70],
		strokecolor=:black, strokewidth=2)
	Label(fig[3, 1][1, 1][2, 1], justification=:left, halign=:left, padding=10, color=:black, word_wrap = true,
		text = "Implikationen:\n" *
			"Die Ausprägung des Suckers’ Payoff hat im iterierten Gefangenendilemma " *
			"eine entscheidene Rolle: Sinkende Werte, insbesondere sobald sie negativ " *
			"werden, schränken die Entwicklung kooperativer Strategien ein oder verhindern " *
			"sie gleich ganz. Es lohnt sich nicht, von defektierenden Strategien abzuweichen, " *
			"da die Strafe einfach zu hoch ist."
	)
end

#-----------------------------------------------------------------------------------------
"""
	plot_history!(fig, sm::SimulationModel)

Create axis for the average strategy history of a single SimulationModel.
"""
function plot_history!(fig, sm::SimulationModel)

	axis = Axis(fig,
		xgridvisible = true, ygridvisible = true,
		xticks = 0:0.2:1, yticks = 0:0.2:1,
		xticklabelsize = 12,
		yticklabelsize = 12,
		aspect = DataAspect()
	)

	lines!( axis, getindex.(sm.history, 1), getindex.(sm.history, 2),
		color = LinRange(1, 10, length(sm.history)), colormap=Reverse(:copper),
		linewidth = 5
	)
	xlims!(axis, 0, 1)
	ylims!(axis, 0, 1)

	text!(
		axis,
		0.05, 0.95,
		text = "S = $(sm.payoffs[1,2])",
		align = (:left, :top),
		space = :relative,
		fontsize = 14,
		color = :black,
	)
end

#-----------------------------------------------------------------------------------------
"""
	present_plots!(fig, models::Vector{SimulationModel}, nCols = 3)

Plot histories of multiple simulation models and arrange them in a
grid with nCols columns, a title and a description of dimensions.
"""
function present_plots!(fig, models::Vector{SimulationModel}, nCols = 3)
	
	# Generate title and subtitle above all plots
	Label(fig[3, 1][1, 2:4][1, 1], color=:black, font=:bold, tellwidth=false,
		text = "Trajectory of Forgiveness against Reciprocation over time"
	)

	# Distribute plots across columns
	nModels = length(models)
	for i in 1:nModels
		row = div(i-1, nCols) + 1
		col = mod(i-1, nCols) + 1
		plot_history!(fig[3, 1][1, 2:4][2, 1][row, col], models[i])
	end

	Label(fig[3, 1][1, 2:4][3, 1], color=:black, tellwidth = false, padding = (0, 0, 0, -15), 
		text = "[x] = Reciprocation, [y] = Forgiveness"
	)
end

#-----------------------------------------------------------------------------------------
# TODO times new roman, sollen wir es verwenden? Können wir die References noch kürzen?
"""
	present_references!( fig)

In the given figure, list all sources of information cited in all the above descriptions.
The journal 'Constructivist Foundations' requires references to be in Times New Roman font.
"""
function present_references!( fig)
	Box(fig[4, 1], color = Makie.cgrad(:cool, alpha = 0.03)[100],
		strokecolor=:black, strokewidth=2, tellwidth = false)
	Label(fig[4, 1], justification=:left, halign=:left, padding=10, color=:black,
		word_wrap = true, fontsize = 11, tellwidth = false,
		text="References: " *
			"Laland, K. (2024) DIY evolution. New Scientist 264(3520): 26–29 / " *
			"Nowak, M.A (2006) Evolutionary dynamics. Harvard University Press / " *
			"Wilson D.S. & Wilson, E.O. (2007)Evolution: Survival of the selfless. " *
			"New Scientist 196(2628): 42–46."
	)
end

#-----------------------------------------------------------------------------------------
# Demo methods:
#-----------------------------------------------------------------------------------------
"""
	demo()

Use PD dynamics to investigate the following research question:
	Can cooperative strategies still prevail, even when the price for cooperation
	in the face of defection increases?
"""
function demo()
	# Arguments for SimulationModel constructor
	payoffs = [4 0; 5 1]
	nStrategies = 100			# 100
	nTrials = 10				# 50
	nMutations = 2000			# 1000
	nGenerations = 150			# 20
	mu = 0.01					# 0.01
	initStrat = [0.075, 0.075]
	initEpsilon = 0.025

	# Suckers' Payoff variation
	sRange = 0.7:-0.2:-0.7

	sVect = collect(sRange)
	nRuns = length(sRange)
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
		setS!(models[run], sVect[run])
		run!(models[run])
	end

	present(models)
end

function demo2()
	# Arguments for SimulationModel constructor
	payoffs = [4 0; 5 1]
	nStrategies = 5				# 100
	nTrials = 5					# 50
	nMutations = 100			# 1000
	nGenerations = 70			# 20
	mu = 0.01					# 0.01
	initStrat = [0.075, 0.075]
	initEpsilon = 0.025

	sm = SimulationModel(
			payoffs,
			nStrategies,
			nTrials,
			nMutations,
			nGenerations,
			mu,
			initStrat,
			initEpsilon
			)

	set!(sm)
	setS!(sm, -100)
	novak!(sm)
	
	display(sm.novak)
	Main.Infiltrator.@infiltrate

end

end		# ... of module SuckersPayoff