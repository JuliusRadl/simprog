#========================================================================================#
"""
	Cooperation

Module Cooperation: A model of players interacting, using strategies of varying cooperation.

This program investigates the following research question:
	By what route can altruistic cooperation infiltrate a thoroughly exploitative population?

Author: Sophia Prögler, Julius Radl, 2026-01-18
"""
module Cooperation

using Random, GLMakie

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Players

A simulation of evolutionary game strategies for PD interactions in a population of nStrategies.
"""
mutable struct Players
	payoff::Matrix{Float64}			# Payoff matrix
	strategies::Vector{Vector}		# Strategies of current population
	x::Vector{Float64}				# Current frequency of strategies
	novak::Matrix{Float64}			# Novak's matrix of expected payoffs

	"""
		Players( payoff, nstrategies)

	The one-and-only Players constructor: Create a Players population of size nstrategies that will
	interact using the given payoff matrix.
	"""
 	function Players( A::Matrix=[4 0;5 1], nstrategies=10)
		new(
			Float64.(A),							# Payoff matrix
			[rand(2) for _ in 1:nstrategies],	    # p and q strategy parameters
			ones(nstrategies)/nstrategies,			# Uniform strategy distribution
			zeros(Float64,2,2),						# Initial novak matrix
		)
	end
end

#-----------------------------------------------------------------------------------------
"""
	Exploration

A structure for saving and modifying parameters for repeated exploration of developmental
trajectories.
"""
mutable struct Exploration
	players::Players
	nTrials::Int64
	nMutations::Int64
	nGenerations::Int64
	mu::Float64
	initStrat::Vector{Float64}
	initEpislon::Float64

	function Exploration(
		players,
		nTrials,
		nMutations,
		nGenerations,
		mu,
		initStrat,
		initEpsilon)

		new(
			players,
			nTrials,
			nMutations,
			nGenerations,
			mu,
			initStrat,
			initEpsilon)
	end
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	mutate!( players, mu)

Replace the least fit fraction mu of the current strategies by new, randomly selected strategies,
and adjust frequencies accordingly.
"""
function mutate!( players::Players, mu::Float64=0.2)
	
    # Legibility
	x = players.x
	nStrategies = length(x)

	# Number of worst Strategies
	nWorstStrats = ceil(Integer, length(x) * mu)

	# Find worst strategies
	ranking = sortperm(x)
	worstStratsIndex = ranking[1:nWorstStrats]
	bestStratsIndex = ranking[nWorstStrats+1:end]

	# Generate new strategies
	players.strategies[worstStratsIndex] = [rand(Float64, 2) for strategy in 1:nWorstStrats]

    # Adjust frequencies
	players.x[worstStratsIndex] .= (1/nStrategies)
	remaining_freq = 1 - sum(players.x[worstStratsIndex])
	players.x[bestStratsIndex] .= remaining_freq * (
		players.x[bestStratsIndex] / sum(players.x[bestStratsIndex]))

	# Recalculate Payoff matrix
	novak!(players)
end

#-----------------------------------------------------------------------------------------
# TODO: Löschen
"""
	vary!(players)

Slightly vary strategies.
"""
function vary!(players::Players, epsilon=0.01)
	strategies = [strat + epsilon*(2*rand(2).-1) for strat in players.strategies]
	players.strategies = map(strategies) do strat
		# Ensure strategy lies within the interval (0,1):
		max.(epsilon,min.(1-epsilon,strat))
	end
end

#-----------------------------------------------------------------------------------------
# TODO: Löschen
"""
	redistribute!(players)

Slightly redistribute population size.
"""
function redistribute!(players::Players, epsilon=0.01)
	nStrategies = length(players.x)
	players.x = players.x + epsilon * (2*rand(nStrategies).-1)
	players.x = players.x ./ sum(players.x)
end

#-----------------------------------------------------------------------------------------
"""
	set!( players, strategy)

Set all current strategies to the given strategy with random variation epsilon.
"""
function set!( players::Players, strategy::Vector{Float64}, epsilon=0.01)
	nstrategies = length(players.x)

	strategies = [copy(strategy) + epsilon*(2*rand(2).-1) for _ in 1:nstrategies]
	players.strategies = map(strategies) do strat
		# Ensure strategy lies within the interval (0,1):
		max.(epsilon,min.(1-epsilon,strat))
	end
	players.x .= 1/nstrategies
end

#-----------------------------------------------------------------------------------------
"""
	simulate!( players, ngenerations)

Simulate the Players dynamics for ngenerations, starting from initial state x0.
"""
function simulate!( players::Players, ngenerations::Int=1)
	x = players.x       # initial state x0
    A = players.novak   # long term payoffs
    dt = 1.0			# step length

    # Calculate population trajectory:
	for _ in ngenerations

		# Perform RK2 half-step:
		R = x' * A * x
		dxdt = x .* (A * x .- R)
		half_step = x + x .* dxdt * (dt/2)
		
		# Perform RK2 full-step:
		new_R = half_step' * A * half_step
		new_dxdt = half_step .* (A * half_step .- new_R)
		x = x + x .* new_dxdt * dt

		# Normalize frequencies
		x = x ./ sum(x)
	end

    players.x = x
end

#-----------------------------------------------------------------------------------------
"""
	novak!( players)

Recalculate Novak's expected payoff from current strategies and payoff matrix.
"""
function novak!( players::Players)

	# Legibility
	A_pd = players.payoff
	S = players.strategies

	# Calculating long term payoffs
	Te = A_pd[2, 1]	# Temptation to defect
	Re = A_pd[1, 1]	# Reward for mutual cooperation
	Pu = A_pd[2, 2]	# Punishment for mutual defection
	Su = A_pd[1, 2]	# Sucker's payoff
	
	p = [pq[1] for pq in S]	# Reciprocation
	q = [pq[2] for pq in S]	# Forgiveness
	r = p .- q

	s_ij = (r .* q' .+ q) ./ (1 .- (r * r'))	# Prob. that S[i] cooperates with S[j] long-term
	s_ji = s_ij'								# The reverse

	Te_part = Te * ((1 .- s_ij) .* s_ji)
	Re_part = Re * (s_ij .* s_ji)
	Pu_part = Pu * ((1 .- s_ij) .* (1 .- s_ji))
	Su_part = Su * (s_ij .* (1 .- s_ji))

	A = Te_part + Re_part + Pu_part + Su_part
	
	# Assigning the long term payoffs
	players.novak = A
end

#-----------------------------------------------------------------------------------------
"""
	avgstrategy( players)

Return average strategy across current population.
"""
function avgstrategy( players::Players)
	players.x'*players.strategies
end

#-----------------------------------------------------------------------------------------
# Display methods:
#-----------------------------------------------------------------------------------------
"""
	behaviour( p::Float64, q::Float64)

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
function behaviour( strategy::Vector)
	p,q = strategy
	
	if q > 1-cutoff
		return (p<cutoff) ? "Ctry" : ((p<1-cutoff) ? "Plac" : "AllC")
	elseif q > cutoff
		return (p<cutoff) ? "Xplt" : ((p<1-cutoff) ? "Gmbl" : "GTFT")
	else
		return (p<cutoff) ? "AllD" : ((p<1-cutoff) ? "STFT" : "RTFT")
	end
end
const cutoff = 0.2						# Strategy classification cutoff

#-----------------------------------------------------------------------------------------
"""
	show( players, bpause)

Display the current Players frequencies at the console in descending frequency order.
"""
function show( players::Players, proportion=1, bpause::Bool=false)
	proportion = max(0,min(1,proportion))
	ntoshow = ceil(Int,proportion*length(players.x))
	toshow = sortperm(players.x,rev=true)[1:ntoshow]
	for i in 1:ntoshow
		println("$(behaviour(players.strategies[toshow[i]])): $(players.x[toshow[i]])")
	end
	println()
	if bpause
		# Pause to allow user to digest output
		readline()
	end
end

#-----------------------------------------------------------------------------------------
"""
	investigate(players)

Explore which strategies survive.
"""
function investigate(Exploration)

	avghistory = [[0.0,0.0] for _ in 1:nMutations+1]	# History averaged over ntrials.

	set!(exp.players, exp.initStrat, exp.initEpsilon)				# Set up initial exploitative population

	for _ in 1:exp.nTrials
		# Perform a single trial simulation:
		set!(exp.players, exp.initStrat, exp.initEpsilon)			# Reset initial exploitative population
		avghistory[1] += avgstrategy(exp.players)

		for mut in 1:exp.nMutations
			# Mutate then simulate nGenerations:
			mutate!(exp.players, exp.mu)						# Replace current weakest fraction of strategies
			simulate!(exp.players, exp.nGenerations)			# Iterate PD interactions
			avghistory[mut+1] += avgstrategy(exp.players)		# Accumulate new simulation results
		end
	end

	avghistory = avghistory ./ nTrials					# Normalize results
end

#-----------------------------------------------------------------------------------------
"""
	varysuckers()


#-----------------------------------------------------------------------------------------
"""
	present(history::Vector)

Present the results of a Players run.
"""
function present(history::Vector)
	fig = Figure(fontsize=16, regular="Helvetica", linewidth=5, size=(1412,1000))

	present_title!(fig)
	present_context!(fig)
	present_problem!(fig)
	present_method!(fig)
	present_results!(fig)
	present_graph!(fig, history)
	present_implications!(fig)
	present_references!(fig)

	fig
	# TODO: später reinkommentieren
	# save("CooperationInfoSheet.pdf",fig)
	# nothing
end

#-----------------------------------------------------------------------------------------
"""
	present_title!( fig)

In the given figure, describe the relevance of the Cooperation model in the context of current
scientific debate.
"""
function present_title!( fig)
	Label(fig[1, 1:2], fontsize = 30, justification=:left, halign=:left, color=:darkgreen, font=:bold,
		text = "Is group selection essential to the evolution of cooperation?"
	)
	Label(fig[2, 1:2], fontsize = 20, justification=:left, halign=:left, color=:darkgreen, font=:italic,
		text = "A computational game dynamics study in evolutionary biology\n" *
			"Niall Palfreyman, Weihenstephan-Triesdorf University of Applied Sciences, 3.11.2025"
	)
end

#-----------------------------------------------------------------------------------------
"""
	present_context!( fig)

In the given figure, describe the relevance of the Cooperation model in the context of current
scientific debate.
"""
function present_context!( fig)
	Box(fig[1:2,3], cornerradius=10, color=:darkgreen, strokecolor=:darkgreen, strokewidth=2)
	Label(fig[1:2,3], justification=:left, halign=:left, padding=10, color=:white,
		text = "Context:\n" *
			"Wilson & Wilson (2007) argue that group\n" *
			"selection drives evolution of cooperative\n" *
			"behaviour. However, many theorists claim\n" *
			"group selection is irrelevant to evolution."
	)
end

#-----------------------------------------------------------------------------------------
"""
	present_problem!( fig)

In the given figure, describe the specific problem that the Cooperation model addresses in current
scientific debate.
"""
function present_problem!( fig)
	Box(fig[3,1], cornerradius=10, color=:blue, strokecolor=:darkgreen, strokewidth=2)
	Label(fig[3,1], justification=:left, halign=:left, padding=10, color=:white,
		text = "Problem:\n" *
			"We observe cooperative behaviour\n" *
			"in the biological world; however, in\n" *
			"such examples, it is often unclear\n" *
			"whether these behaviours originally\n" *
			"evolved through group selection."
	)
end

#-----------------------------------------------------------------------------------------
"""
	present_method!( fig)

In the given figure, describe the exact, algorithmic sequence of steps and measurements that others
scientists must use in order to reproduce your findings.
"""
function present_method!( fig)
	Box(fig[3,2:3], cornerradius=10, color=:blue, strokecolor=:darkgreen, strokewidth=2)
	Label(fig[3,2:3], justification=:left, halign=:left, padding=10, color=:white,
		text = "Method:\n" *
			"Group selection presupposes that the subjects of evolution are niche-" *
			"constructing (Laland 2024) developmental (Puentedura 2007)\n" *
			"processes (evo-eco-devo). Here, we simulate computationally (Nowak 2006)" *
			"a population of Prisoner’s Dilemma (PD) processes whose\n" *
			"strategy can mutate arbitrarily with low probability. Strategies are selectively " *
			"punished for exhibiting (with probability q) the cooperative\n" *
			"behaviour of forgiveness, unless this forgiveness is reciprocated " *
			"(with probability p) by other players. In this case, selection rewards the\n" *
			"mutual interaction between forgiveness and reciprocation.\n" *
			"We discuss whether this dynamical system demonstrates that group selection is " *
			"essential to the evolution of cooperation."
	)
end


#-----------------------------------------------------------------------------------------
"""
	present_results!( fig)

In the given figure, describe the precise, uninterpreted and error-prone empirical data that
arose from implementing the method.
"""
function present_results!( fig)
	Box(fig[4,1], cornerradius=10, color=:red, strokecolor=:darkgreen, strokewidth=2)
	Label(fig[4,1], justification=:left, halign=:left, padding=10, color=:white,
		text = "Results:\n" *
			"In the graph, cooperation evolves in\n" *
			"a mutating population of PD defectors.\n" *
			"Initially, forgiving behaviour is\n" *
			"penalised: its representation in the\n" *
			"population falls as reciprocating rises.\n" *
			"Subsequently, forgiving rises to q≈0.6\n" *
			"and reciprocating rises towards p≈1.0."
	)
end

#-----------------------------------------------------------------------------------------

function present_graph!(fig, history)
	
	axis = Axis(fig[4:5, 2][2, 2],
		title = "Trajectory of Forgiveness against Reciprocation over time",
		xlabel = "Reciprocation", ylabel = "Forgiveness",
		xgridvisible = true, ygridvisible = true,
		xticks = 0:0.1:1, yticks = 0:0.1:1,
		backgroundcolor=(:green, 0.1)
	)

	scatterlines!( axis, getindex.(history,1), getindex.(history,2),
		color=:blue,
		markercolor=:red, markersize=10
	)
	xlims!(axis, 0, 1)
	ylims!(axis, 0, 1)
end

#-----------------------------------------------------------------------------------------

function present_multiple_graphs!(fig, history)
	
	axis = Axis(fig[4:5, 2][2, 2],
		title = "Trajectory of Forgiveness against Reciprocation over time",
		xlabel = "Reciprocation", ylabel = "Forgiveness",
		xgridvisible = true, ygridvisible = true,
		xticks = 0:0.1:1, yticks = 0:0.1:1,
		backgroundcolor=(:green, 0.1)
	)

	scatterlines!( axis, getindex.(history,1), getindex.(history,2),
		color=:blue,
		markercolor=:red, markersize=10
	)
	xlims!(axis, 0, 1)
	ylims!(axis, 0, 1)
end

#-----------------------------------------------------------------------------------------
"""
	present_implications!( fig)

In the given figure, describe the meaning and implications of your results in relation to the
previous description of the problem and its context.
"""
function present_implications!( fig)
	Box(fig[5,1], cornerradius=10, color=:red, strokecolor=:darkgreen, strokewidth=2)
	Label(fig[5,1], justification=:left, halign=:left, padding=10, color=:white,
		text = "Implications:\n" *
			"Cooperation evolves as defectors\n" *
			"construct around themselves a\n" *
			"mutually reciprocating niche of\n" *
			"players that reward forgiveness.\n" *
			"Since this niche is selected due to\n" *
			"its group property of reciprocating\n" *
			"behaviour between its members,\n" *
			"our results suggest that selection of\n" *
			"groups is important for the evolution\n" *
			"of cooperation."
	)
end

#-----------------------------------------------------------------------------------------
"""
	present_references!( fig)

In the given figure, list all sources of information cited all the above descriptions. The journal
Constructivist Foundations requires references to be in Times New Roman font.
"""
function present_references!( fig)
	Box(fig[4:5,3], cornerradius=10, color=(:lime, 0.5), strokecolor=:darkgreen, strokewidth=2)
	Label(fig[4:5,3], font="Times New Roman", fontsize=20, justification=:left, halign=:left,
		padding=10, color=:black,
		text="References:\n" *
			"•   Laland, K. (2024) DIY evolution.\n    New Scientist 264(3520): 26–29.\n" *
			"•   Nowak, M.A (2006) Evolutionary\n    dynamics. Harvard University\n    Press.\n" *
			"•   Puentedura, R.R. (2007) The\n    Baldwin effect in the age of\n    computation. " *
			"In: Weber, B.H. &\n    Depew, D.J. (eds): Evolution and\n    learning. MIT Press.\n" *
			"•   Wilson D.S. & Wilson, E.O. (2007)\n    Evolution: Survival of the selfless.\n" *
			"    New Scientist 196(2628): 42–46."
	)
end

#-----------------------------------------------------------------------------------------
# Demo methods:
#-----------------------------------------------------------------------------------------
"""
	demo()

Use PD dynamics to investigate the following research question:
By what route can altruism infiltrate a thoroughly exploitative population?
"""
function demo()
	payoffs = [4 0; 5 1]
	nStrategies = 100
	nTrials = 50
	nMutations = 100
	nGenerations = 20
	mu = 0.01
	initStrat = [0.075, 0.075]
	initEpsilon = 0.025

	history = investigate(
		payoffs,
		nStrategies,
		nTrials,
		nMutations,
		nGenerations,
		mu,
		initStrat,
		initEpsilon
	)

	present(history)
end

end		# ... of module Players