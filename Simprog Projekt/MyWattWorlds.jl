#========================================================================================#
"""
	WattWorlds

Module WattWorlds: A model of narrative stabilisation to demonstrate narrative stabilisation. The
WattWorld structure comprises a homogenous population of individually narrative Players, each
characterised by a transient state a (0 < A_LWB < a < A_UPB) and its behaviour parameters B. Both
behaviour parameters are constrained to lie between bounds:
	B ∈ (0 < B_LWB < abs(B[1]) < B_UPB) × (0 < C_LWB < B[2] < C_UPB).

Players' actions are coordinated by their common access to a single, globally available resource R.
Each Player influences the value of R in proportion to its (positive-valued) activation state a,
either producing or consuming R according to the value of its behaviour parameters B. B specifies
Hill-function parameters:
	-	B is a Michaelis-Menten-style half-saturation parameter specifying the threshold value of R
		beyond which the Player's consumption/production rate rapidly rises from 0 to 1;
	-	C is the corresponding cooperation parameter, specifying the abruptness of this threshold.
	
Each Player functions as an integral controller of R that, depending upon its B-value, potentially
contributes one rein of an integral-rein Watt governor.

Structural mutation of WattWorld occurs through stochastic variation of Players' B-values, the
extent of this variation depending rheolectically upon the current coordination value R. In the
simple WattWorld model, narrative stabilisation is trivial: B-variability is zero when R=1.0, rising
monotonically to a maximum value of 1.0 as R moves further away from the value 1.0.

The WattWorld model explores two research questions relating to this simple rheolectic function:
	a) Can narrative stabilisation lead to stabilisation of the WattWorld narrative, and if so,
	b) Does this stabilised narrative enact the dynamics of engagement with an exogenously
		determined, time-dependent injection supply(t) of the resource R?

Author: Niall Palfreyman, 15/09/2024
Revised: 15/09/2025
"""
module WattWorlds

include("./AgentTools.jl")

using GLMakie, Random, Agents, Statistics, .AgentTools

#-----------------------------------------------------------------------------------------
# Module constants:
#-----------------------------------------------------------------------------------------
# TODO string above constant is displayed in tooltip for constant

"Radius in which to look for agents when calculating inhibition per cell"
const NBR_RADIUS = 2
"Extend of the continuous ABM Space"
const EXTENT = (50, 50)
"Duration of simulation"
const DURATION = 15e6
"Step length for RK2 simulation"
const RK2_STEP = 2.0
"Stride length for graphical data display compression"
const GRAPHICS_COMPRESSION = 50000
"Maximum magnitude of stochastic B-variation"
const DELTA = 0.05
"Lower bound of Players' activation state a"
const A_LWB = 1e-300
"Upper bound of Players' activation state a"
const A_UPB = 3.0
"Minimum level of activation state a for Players to be graphically reported"
const A_REPORT = 0.01
"Lower bound of Players' Hill-saturation parameter abs(B.K)"
const K_LWB = 0.1
"Upper bound of Players' Hill-saturation parameter abs(B.K)"
const K_UPB = 2.0
"Lower bound of Players' Hill-cooperativity parameter B.n"
const n_LWB = 1.0
"Upper bound of Players' Hill-cooperativity parameter B.n"
const n_UPB = 9.0
"Size of domain of Players' Hill-cooperativity parameter B.n"
const n_SPAN = n_UPB - n_LWB
"Upper bound of WattWorld resource value R"
const R_UPB = 10.0
"Depletion constant for Players' activation state a"
const BETA_A = 0.1
"Depletion constant for resource R"
const BETA_R = 1.0
"Radius of stability well"
const CAPTURE_RADIUS = 1.0
"Monomial steepness of stability well"
const CAPTURE_CONCAVITY = 2.5
"Descriptions of available benchmark test regimes"
const REGIMES = [
	"Basic stabilisation (constant feed and context)"
	"Ontogenic adjustment (single-step feed; constant context)"
	"Phylogenic adjustment (constant feed; single-step context)"
	"Ontogenic stabilisation (periodic feed; constant context)"
	"Phylogenic stabilisation (constant feed; periodic context)"
	"Agency (Transfer fast periodic feed control across slow periodic context)"
]

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Player

The Player type encapsulates a single WattWorld player, each characterised by its:

	* Activation state: a ∈ [A_LWB,A_UPB]
	* Behavioural half-saturation constant B.K ∈ [-B_UPB,-B_LWB] U [B_LWB,B_UPB]
	* Behavioural cooperativity constant B.n ∈ [0,n_UPB]
"""
# Agent in 2-dimensional continuous space (automatically mutable using @agent)
@agent struct Player(ContinuousAgent{2, Float64})
	a::Float64										# Activation state variable
	B::NamedTuple{(:K,:n), Tuple{Float64,Float64}}	# Behavioural Hill constants
	groupcolor::Symbol								# color indicating producer or consumer
end

#-----------------------------------------------------------------------------------------
# """
# 	WattWorld

# A WattWorld is a collection of N Players playing a shared, iterated integral-rein game.
# """
# umwandeln in abm-passende initialiserungsfunktion
# 
# mutable struct WattWorld
# 	players::Vector{Player}						# Vector of Players
# 	R::Float64									# Global resource level accessible to Players
# 	feed::Float64								# Current feed rate of resource R
# 	omega::Float64								# Target resource level at which ΔB=0
# 	regime::Int									# Feed/omega regime for running WattWorld
# 	t::Float64									# Current time

# 	# Private backup fields for implementation of Runge-Kutta-2 step!():
# 	backup_a::Vector{Float64}
# 	backup_R::Float64
# 	backup_t::Float64

# 	# Unique constructor:
# 	function WattWorld( n_players=N_PLAYERS; regime=1)


# 		new(
# 			[Player() for _ in 1:n_players],	# N random Players
# 			A_LWB,								# Arbitrarily tiny initial value of R
# 			0.2,								# Resource supply rate
# 			1.0,								# Omega stabilisation fixed point of R
# 			regime,								# Use this feed/omega regime
# 			0.0,								# Current time set to zero
# 			zeros(Float64,n_players),			# N placeholders for backup activations
# 			0.0,								# Placeholder for backup resource value
# 			0.0									# Placeholder for backup time value
# 		)
# 	end
# end

#-----------------------------------------------------------------------------------------
function wattWorld(;		# keyword arguments follow, setting defaults
							# these arguments can be changed with sliders
							# so not everything goes in here
	feed	= 0.2,			# Current feed rate of resource R
	omega	= 1.0,			# Target resource level at which ΔB=0
	regime	= 1,			# Feed/omega regime for running WattWorld
	n_players = 5			# Number of Players
)
	# house keeping
	@assert n_players > 1					# Number of Players at least 2
	@assert 1 <= regime <= length(REGIMES)	# Valid regime number
	GC.gc()									# Clear up garbage from previous run
	Random.seed!(5)							# Make simulation reproducible

	# set properties as a dictionary, later to be called with model.<property>
	properties = Dict(
		# fields that can be changed with sliders
		:feed => feed,
		:omega => omega,
		:regime => regime,
		:n_players => n_players,

		# initialise space
		:extent => EXTENT,
		:R => fill(A_LWB, EXTENT),		# Initialize resource matrix with tiny starting value

		# initial number of steps taken is 0
		:t => 0,
		:dt => RK2_STEP,

		# backup fields for RK2-Step
		:backup_a => zeros(Float64, n_players),	# N placeholders for backup activations
		:backup_R => 0.0,						# Placeholder for backup resource value
		:backup_t => 0						# Placeholder for backup time value
	)

	# create model with Agents function StandardABM
	model = StandardABM(
		Player, ContinuousSpace(EXTENT; spacing=1.0);
		properties, agent_step!, model_step!
	)

	# Fill model with Players
	for _ in 1:n_players
		K = wrap((2rand()-1)*K_UPB)		# Arbitrary valid half-saturation ...
		n = n_LWB + rand()*n_SPAN		# and cooperativity
		theta = 2π*rand()				# Random angle
		player = add_agent!(model, (cos(theta), sin(theta)), # random angle
			A_LWB, # Initial tiny activation
			(K = K, n = n), # structure variables / genome
			:white	# temporary group identity color
		)
		set_groupcolor!(player)
	end

	# return the model
	return model
end

#-----------------------------------------------------------------------------------------
"""
	Snapshot

A Snapshot records the essential state of a WattWorld for later trajectory plotting.
"""
struct Snapshot
	players::Vector{Player}						# WattWorld's Players
	R::Float64									# WattWorld's global resource level
	F::Float64									# WattWorld's current feed rate
	Ω::Float64									# WattWorld's current omega point
	t::Float64									# WattWorld's current time

	# Unique copy constructor
	function Snapshot(model)
		new( deepcopy(allagents(model)), model.R, model.feed, model.omega, model.t)
	end
end

#-----------------------------------------------------------------------------------------
# Agent related methods
#-----------------------------------------------------------------------------------------
"""
	agent_step!(player, model)

"""
function agent_step!(player, model)
	# check resource level at location

	# move

	# decrease own energy level

	# try to reproduce

	# increasing or decreasing the resource is part of model_step!

	return
end

#-----------------------------------------------------------------------------------------
"""
	set_groupcolor!(player)

Set players color to indicate whether they consume or produce resource.
"""
function set_groupcolor!(player)
	# TODO make more fine grain colors for different K values
	if player.B.K > 0
		player.groupcolor = :blue # consumers with positive k's lower ressource => make them blue
	else
		player.groupcolor = :red
	end
end

#-----------------------------------------------------------------------------------------
# Model related methods
#-----------------------------------------------------------------------------------------
"""
	model_step!(model)


"""
function model_step!(model)
	# set omega and feed according to current regime
	apply_regime!(model)

	# update a's and R's
	rktwo_step!(model)

	# diffuse R's

	
end

#-----------------------------------------------------------------------------------------
"""
	rk_backup!( model)

Backup WattWorld's activation, resource and steps taken into its private backup fields for RK2 halfstep.
"""
function rk_backup!( model)
	model.backup_a[:] = (p->p.a).(allagents(model))
	model.backup_R    = model.R
	model.backup_t    = model.t
end

#-----------------------------------------------------------------------------------------
"""
	rk_restore!( model)

Restore WattWorld's activation, resource and steps taken from its private backup fields for RK2 halfstep.
"""
function rk_restore!( model)
	# pairs automatically creates an index, even if the argument is not a dict
	for (i,p) in pairs(allagents(model))
		p.a = model.backup_a[i]
	end
	model.R = model.backup_R
	model.t = model.backup_t
end

#-----------------------------------------------------------------------------------------
"""

Calculates a vector of activation RoC for all players.
"""
function calculate_rocs(model)
	acts = (p->p.a).(allagents(model))						# Activations of Players
	Ks = (p->p.B.K).(allagents(model))						# Half-saturation constants K of Players

	# differential equation for RoC of a: inhibition and activation terms
	p_idxs = get_player_idxs(model)

	# Get vector containing inhibition for each player, based on the cell they are in
	# (precalculate inhibs for each cell)
	cell_idxs = get_player_cell_idxs(p_idx)

	# --------------- da ------------------
	inhib_mat = get_cell_inhibitions(model, cell_idxs)
	inhibs = get_player_inhibs(inhib_mat, p_idxs) # Players' local action inhibits own activation

	# Get vector of R values for each player
	R_vect = R[[p_idxs]]

	# Get vector containing individual activation for each player
	activation = ((R, K)->saturation(R,K)).(R_vect, Ks)			# R-dependent activation of Players
	
	# net activation
	net_activation = inhibs .* activation			# Net growth of each Player's activation

	# full differential equations
	da = (net_activation .- BETA_A) .* acts				# Growth - depletion of activations

	# ----------- dR --------------------
	# initialise resource use matrix (sum of a * K per per player in cell)
	res_use = zeros(model.extent)

	# calculate resource use for cells containing players
	for (p, idx) in pairs(allagents(model), p_idxs)
		res_use[idx] += (p.a * p.B.K)
	end
	
	# calculate RoC matrix for Resource R
	dR = model.feed .- BETA_R .* model.R .- res_use	# Feed rate minus depletion and consumption/production

	return (da, dR)
end

#-----------------------------------------------------------------------------------------
"""

Get vector containing indices of each player.
"""
function get_player_idxs(model)
	# allocate empty vector of tuples
	agents = allagents(model)
	p_idxs = Vector{Tuple{Int,Int}}(undef, length(agents))
	for (i, p) in pairs(agents)
		# needs reference matrix for size
		# and model to interpret geometry (like continous space)
		idx = get_spatial_index(p.pos, model.extent, model)
		p_idxs[i] = idx
	end
	return p_idxs
end

#-----------------------------------------------------------------------------------------
"""

Get vector of inhibitions for each player from an inhibition matrix.
"""
function get_player_inhibs(inhib_mat, player_idxs)
	p_inhibs = inhib_mat[player_idxs]
	return p_inhibs
end

#-----------------------------------------------------------------------------------------
"""

Gets indices of all cells that contain atleast one player.
"""
function get_player_cell_idxs(player_idxs)
	p_cell_idxs = unique(player_idxs)
	return p_cell_idxs
end

#-----------------------------------------------------------------------------------------
"""

Calculates inhibition per cell based on the activation and saturation parameters of
agents in neighbouring cells. Avoids calculating inhibition per agent (and checking each agents
neighbours), but loses some accuracy. Should help performance for clumping agents.
accuracy.

player_cell_idx must be vector of integer indices that fit model extent.
"""
function get_cell_inhibitions(model, player_cell_idxs)
	# initialise inhibitions
	inhib_mat = ones(model.extent) # no inhibition for empty cells
	
	# iterate over all cells that contain players
	for idx in player_cell_idxs
		# get all agents in this cell and adjacent cells
		# respects periodic boundaries
		cell_center = (idx[1] - 0.5, idx[2] - 0.5)
		# TODO could be faster instead, check later:
		# x = sum(nbr.a * nbr.B.K for nbr in nearby_agents(cell_center, model, NBR_RADIUS))
		nbrs_iter = nearby_agents(cell_center, model, NBR_RADIUS)
		nbrs = collect(nbrs_iter)
		
		# calculate inhibition based on neighbours
		acts = (nbr->nbr.a).(nbrs)
		Ks = (nbr->nbr.B.K).(nbrs)
		
		# hill saturation function
		x = sum(acts .* Ks)
		inhib = saturation(x, -1.0)
		
		# set in inhibition matrix
		inhib_mat[idx] = inhib
	end

	return inhib_mat
end

#-----------------------------------------------------------------------------------------
"""
	rktwo_step!( model)

Develop WattWorld through a single Runge-Kutta-2 time-step according to model properties.
"""
function rktwo_step!( model)
	# RK2 half-step:
	# backup current parameters
	rk_backup!(model)

	# calculate half step length
	dt_half = model.dt/2.0

	# calculate RoC for a and R
	da,dR = calculate_rocs(model)

	# Add halved RoCs and time step to wattworld parameters
	# use for loop for writing to each player (we are writing into the model,
	# not into temporary variables
	for (i,p) in pairs(allagents(model))
		p.a += da[i] * dt_half
	end
	model.R .+= dR .* dt_half
	model.t .+= dt_half

	# RK2 full-step:
	# use new parameters to caculate new RoCs
	da,dR = calculate_rocs(model)

	# restore backed up parameters for full step
	rk_restore!(model)

	# perform full step with full step length
	for (i,p) in pairs(allagents(model))
		p.a += da[i] * model.dt
	end
	model.R .+= dR * model.dt
	model.t .+= round(model.dt,digits=3)
end

#-----------------------------------------------------------------------------------------
"""
	vary!( model)

With probability given by WattWorld's current instability, vary the behavioural parameters K and n
of all WattWorld Players by up to an amount DELTA within the respective validity domain.
"""
function vary!( model)
	instability = 1 - stability(model)
	# generate deltas as a 2xlength(players) matrix
	deltas = DELTA * instability * (2rand(2,length(allagents(model))).-1)
	for (i,p) in pairs(allagents(model))
		p.B = (
			# limit change of K's within K Domain
			K = wrap(p.B.K + deltas[1,i]),
			# limit change of n's
			n = n_LWB + rem( n_SPAN + p.B.n - n_LWB + deltas[2,i], n_SPAN)
		)
	end
end

#-----------------------------------------------------------------------------------------
"""
	stability( model)

Determine the stability of WattWorld determined by the current resource level R. Currently,
stability is the closeness of R to WattWorld's omega value within a capture well of width
CAPTURE_RADIUS and with monomial degree CAPTURE_CONCAVITY.
"""
function stability( model)
	1.0 - min( 1.0, (abs(model.R-model.omega)/CAPTURE_RADIUS)^CAPTURE_CONCAVITY)
end

#-----------------------------------------------------------------------------------------
"""
	apply_regime!( model)

Apply WattWorld's current feed/omega test regime for benchmarking.
"""
function apply_regime!( model)
	if model.regime <= 1
		# Constant feed and omega:
		model.omega = 1.0
		model.feed = 0.2
	elseif model.regime == 2
		# Single step feed:
		model.omega = 1.0
		model.feed = iseven(2model.t÷DURATION) ? 0.2 : 1.75
	elseif model.regime == 3
		# Single step omega:
		model.omega = iseven(2model.t÷DURATION) ? 1.0 : 2.0
		model.feed = 0.2
	elseif model.regime == 4
		# Slow periodic feed:
		model.omega = 1.0
		model.feed = iseven(model.t÷1e6) ? 0.2 : 1.75
	elseif model.regime == 5
		# Slow periodic omega:
		model.omega = iseven(model.t÷1e6) ? 1.0 : 2.0
		model.feed = 0.2
	else
		# Fast periodic feed with slow periodic omega:
		model.omega = iseven(model.t÷8e5) ? 1.0 : 2.0
		model.feed = iseven(model.t÷2e5) ? 0.2 : 1.75
	end
end

#-----------------------------------------------------------------------------------------
"""
	regime_string( model)

Return a short string description of WattWorld's current regime.
"""
function regime_string( model)
	REGIMES[model.regime]
end

#-----------------------------------------------------------------------------------------
# Simulation / Plotting related methods
#-----------------------------------------------------------------------------------------
"""
	trajectory( model, T; Δt=T/500) → Snapshot

Develop WattWorld forward from the current instant through given time duration T.
"""
function trajectory( model, T::Real; Δt=1.0)

	# get number of time steps
	num_ts = Int(round(T/Δt))+1									# Allow for roundng errors in t
	
	# which time do we arrive at after developing for duration T?
	tfinal = round(model.t+T,digits=3)

	# creates step range between starting time and end time (after duration T)
	ts = range(model.t,tfinal,length=num_ts)
	# creates delta t as step length of that range
	Δt = step(ts)

	# initializes snapshots, one for each time step in the step range
	snap = Snapshot( model)										# Initialise snapshot storage
	snaps = similar([snap],num_ts)
	snaps[1] = snap


	for i in 2:num_ts
		vary!(model)												# Vary watt world if it is unstable
		step!(model,Δt)											# Take watt forwards 1 step

		# Repair any invalid values:

		# bound R to [0, 9]
		model.R = max(0.0,min(R_UPB,model.R))
		
		# replace all players with a's outside [0, 3]
		for (j,p) in pairs(allagents(model))
			if p.a < 0 || p.a > A_UPB
				allagents(model)[j] = Player()
			end
		end
	
		snaps[i] = Snapshot( model)								# Store new Snapshot
	end

	return snaps
end

#-----------------------------------------------------------------------------------------
"""
	report( model)

Report current status of the WattWorld.
"""
# TODO adaptieren oder löschen
function report( model)
	println(
		"Time=$(round(model.t,digits=3)), Resource=$(model.R), Action=$(action(model)), ",
		"Stability=$(stability(model)):"
	)
	for (i,p) in pairs(allagents(model))
		if p.a > A_REPORT
			println( "    Player $i: B=$(p.B), a=$(p.a)")
		end
	end
end

#-----------------------------------------------------------------------------------------
"""
    mean_ressource_error(model)

Return the mean absolute difference between local ressource level and omega for plotting.
"""
function mean_temperature_error(model)
    return Statistics.mean(abs.(model.R .- model.omega))
end

#-----------------------------------------------------------------------------------------
# Module-level utilities:
#-----------------------------------------------------------------------------------------
"""
	saturation( s::Real, K=1.0, n=1.0)

Return the Hill saturation value of a signal s, half-saturation constant K and cooperativity n.
"""
function saturation( s, K::Real=1.0, n::Float64=1.0)
	if K == 0
		return 0.5
	end

	abs_K = abs(K)
	if n != 1
		abs_K = abs_K^n
		if s < 0 error("Help! s=$s; K=$K; n=$n") end
		s = s.^n
	end

	abs_saturation = abs_K ./ (abs_K .+ s)
	
	(K < 0) ? abs_saturation : (1 .- abs_saturation)
end

#-----------------------------------------------------------------------------------------
"""
	wrap( b::Real)

	Wrap the number b into the permissible B-range [-K_UPB,-K_LWB] U [K_LWB,K_UPB].
"""
function wrap( b::Real)
	span = float(K_UPB - K_LWB)

	if -K_UPB<=b<=-K_LWB || K_LWB<=b<=K_UPB
		return float(b)
	elseif -K_LWB<b<0.0
		return float(b + 2K_LWB)
	elseif 0.0<b<K_LWB
		return float(b - 2K_LWB)
	elseif K_UPB<b<K_UPB+span
		return float(b - 2K_UPB)
	elseif -(K_UPB+span)<b<-K_UPB
		return float(b + 2K_UPB)
	end

	B = mod((b+span),2span)-span
	(B<0) ? span-K_LWB : span+K_LWB
end

#-----------------------------------------------------------------------------------------
# Demos, Publishing, etc.
#-----------------------------------------------------------------------------------------

function demo(regime::Int=1)
	# set up slider ranges
	params = Dict(
		:feed => 0.2:0.001:1.75,
		:omega => 1:0.1:2.0,
		:regime => 1:1:5,
		:n_players => 10:1:100
	)

	# define vector of (anonymous) plot data functions
	mdata = [
		mean_temperature_error
	]

	# define plot labels
	mlabels = [
		"Mean deviation of R from Omega"
	]

	# define plotting parameters
	plotkwargs = (
		agent_color=(player->player.groupcolor),	# group identidy (producer / consumer)
		agent_size=20,
		agent_marker = :circle,
		heatarray = (model->model.R),	# resource level as background colouring
        heatkwargs = (colormap = :thermal, colorrange = (0.0, R_UPB)), # Distribute color in range of possible R values
		mdata = mdata,
		mlabels = mlabels,
		enable_inspection = false
	)

	# get figure and observable
	playground, = abmplayground(wattWorld; params, plotkwargs...)

	# display interactive figure
	playground
end

#-----------------------------------------------------------------------------------------
# """
# 	demo_old()

# Build and run WattWorld.
# """
# function demo_old( regime::Int=1)
# 	# setup watt world with specific regime
# 	model = WattWorld(N_PLAYERS,regime=regime)

# 	# print header string
# 	println(
# 		"\n======= Simulating a $N_PLAYERS-player WattWorld using regime $regime: ",
# 		"$(regime_string(model)) ======="
# 	)

# 	# print report on initial conditions
# 	report(model)

# 	# run wattworld and save snapshots
# 	snapshots = trajectory(model,DURATION,Δt=RK2_STEP)

# 	# print final report
# 	report(model)

# 	# Display resource, activation and behaviour parameters graphically:
# 	# make a figure
# 	fig = Figure(fontsize=30,linewidth=5)
# 	# select a limited number of snapshots
# 	compressed_t = 1:GRAPHICS_COMPRESSION:length(snapshots)
# 	# select values for time axis as x axis
# 	t_axis = ((s->s.t).(snapshots))[compressed_t]
# 	# create axis
# 	ax_R = Axis(fig[1,1], xlabel="time", title="Regime $regime: Global behaviour")
# 	# create line plots for model parameters from snapshots
# 	lines!( ax_R, t_axis, ((s->s.R).(snapshots))[compressed_t], label="Resource")
# 	lines!( ax_R, t_axis, ((s->s.Ω).(snapshots))[compressed_t], label="Omega")
# 	lines!( ax_R, t_axis, ((s->s.F).(snapshots))[compressed_t], label="Feed rate")
# 	# create legend in separate axis
# 	Legend( fig[1,2], ax_R)

# 	# create axis for player parameters
# 	ax_a = Axis(fig[1,3], xlabel="time", title="Activation (a)")
# 	ax_K = Axis(fig[2,1], xlabel="time", title="Half-saturation (B.K)")
# 	ax_n = Axis(fig[2,3], xlabel="time", title="Cooperativity (B.n)")

# 	for i in 1:N_PLAYERS
# 		# get activity levels 
# 		acts = ((s->s.players[i].a).(snapshots))[compressed_t]
# 		# only select players with significant activation in the last half of simulation
# 		if any( acts[end÷2:end] .> A_REPORT)
# 			# Display players with significant activation in the last 1/4 of the simulation:
# 			# plots into the three axes created earlier
# 			lines!( ax_a, t_axis, acts, label="Player $(i)")
# 			lines!( ax_K, t_axis, ((s->s.players[i].B.K).(snapshots))[compressed_t], 
# 				label="Player $(i)"
# 			)
# 			lines!( ax_n, t_axis, ((s->s.players[i].B.n).(snapshots))[compressed_t],
# 				label="Player $(i)"
# 			)
# 		end
# 	end
# 	Legend( fig[2,2], ax_a)
# 	display( fig)
# 	fig
# end

#-----------------------------------------------------------------------------------------
"""
	publish()

Build and run WattWorld, and save results to the publishable graphics used in
Constructivist Foundations article.
"""
function publish( regime::Int=1)
	model = WattWorld(N_PLAYERS,regime=regime)
	println(
		"\n======= Simulating a $N_PLAYERS-player WattWorld using regime $regime: ",
		"$(regime_string(model)) ======="
	)
	report(model)
	snapshots = trajectory(model,DURATION,Δt=RK2_STEP)
	report(model)

	# Display resource, activation and behaviour parameters graphically:
	fig = Figure( fontsize=30,linewidth=5,size=(1500,1200))
	compressed_t = 1:GRAPHICS_COMPRESSION:length(snapshots)
	t_axis = ((s->s.t).(snapshots))[compressed_t]
	ax_R = Axis(fig[1,1], xlabel="time", title="Regime $regime: R, Ω, F(t)")
	lines!( ax_R, t_axis, ((s->s.R).(snapshots))[compressed_t], label="Resource")
	lines!( ax_R, t_axis, ((s->s.Ω).(snapshots))[compressed_t], label="Omega")
	lines!( ax_R, t_axis, ((s->s.F).(snapshots))[compressed_t], label="Feed rate")
	Legend( fig[1,2], ax_R)

	ax_a = Axis(fig[2,1], xlabel="time", title="Activation (a)")
	ax_K = Axis(fig[3,1], xlabel="time", title="Half-saturation (K)")

	for i in 1:N_PLAYERS
		acts = ((s->s.players[i].a).(snapshots))[compressed_t]
		if any( acts[end÷10:end] .> A_REPORT)
			# Display players with significant activation in the last 2/3 of the simulation:
			lines!( ax_a, t_axis, acts, label="Player $(i)")
			lines!( ax_K, t_axis, ((s->s.players[i].B.K).(snapshots))[compressed_t], 
				label="Player $(i)"
			)
		end
	end
	Legend( fig[2,2], ax_a)
	save( "Figure$regime.pdf", fig)
	nothing
end

#-----------------------------------------------------------------------------------------
"""
	demo_all()

Call demo() for all regimes.
"""
function demo_all()
	for regime in 1:length(REGIMES)
		demo(regime)
	end
end

#-----------------------------------------------------------------------------------------
"""
	publish_all()

Call publish() for all regimes.
"""
function publish_all()
	for regime in 1:length(REGIMES)
		publish(regime)
	end
end

end