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

using CairoMakie, Random

#-----------------------------------------------------------------------------------------
# Module constants:
#-----------------------------------------------------------------------------------------
"Duration of simulation"
const DURATION = 15e6
"Step length for RK2 simulation"
const RK2_STEP = 2.0
"Stride length for graphical data display compression"
const GRAPHICS_COMPRESSION = 50000
"Maximum magnitude of stochastic B-variation"
const DELTA = 0.05
"Default number of WattWorld Players"
const N_PLAYERS = 5
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
mutable struct Player
	a::Float64										# Activation state variable
	B::NamedTuple{(:K,:n), Tuple{Float64,Float64}}	# Behavioural Hill constants

	function Player()								# Unique void constructor
		new(
			A_LWB,									# Tiny valid activation value
			(
				K = wrap((2rand()-1)*K_UPB),		# Arbitrary valid half-saturation ...
				n = n_LWB + rand()*n_SPAN			# and cooperativity
			)
		)
	end
end

#-----------------------------------------------------------------------------------------
"""
	WattWorld

A WattWorld is a collection of N Players playing a shared, iterated integral-rein game.
"""
mutable struct WattWorld
	players::Vector{Player}						# Vector of Players
	R::Float64									# Global resource level accessible to Players
	feed::Float64								# Current feed rate of resource R
	omega::Float64								# Target resource level at which ΔB=0
	regime::Int									# Feed/omega regime for running WattWorld
	t::Float64									# Current time

	# Private backup fields for implementation of Runge-Kutta-2 step!():
	backup_a::Vector{Float64}
	backup_R::Float64
	backup_t::Float64

	# Unique constructor:
	function WattWorld( n_players=N_PLAYERS; regime=1)
		@assert n_players > 1					# Number of Players at least 2
		@assert 1 <= regime <= length(REGIMES)	# Valid regime number
		GC.gc()									# Clear up garbage from previous run
		Random.seed!(5)							# Make simulation reproducible

		new(
			[Player() for _ in 1:n_players],	# N random Players
			A_LWB,								# Arbitrarily tiny initial value of R
			0.2,								# Resource supply rate
			1.0,								# Omega stabilisation fixed point of R
			regime,								# Use this feed/omega regime
			0.0,								# Current time set to zero
			zeros(Float64,n_players),			# N placeholders for backup activations
			0.0,								# Placeholder for backup resource value
			0.0									# Placeholder for backup time value
		)
	end
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

	# Unique copy constructor:
	function Snapshot( watt::WattWorld)
		new( deepcopy( watt.players), watt.R, watt.feed, watt.omega, watt.t)
	end
end

#-----------------------------------------------------------------------------------------
# Player methods:
#-----------------------------------------------------------------------------------------
"""
	action( p::Player)

Calculate the action of this Player. A Player's action is a measure of how much influence the
Player has on the resource level R of its WattWorld.
"""
function action( p::Player)
	abs(p.B.K) * p.a
end

#-----------------------------------------------------------------------------------------
# WattWorld methods:
#-----------------------------------------------------------------------------------------
"""
	action( watt::WattWorld)

Calculate total action of all Players in this WattWorld.
"""
function action( watt::WattWorld)
	sum(action.(watt.players))
end

#-----------------------------------------------------------------------------------------
"""
	rk_backup!( watt::WattWorld)

Backup WattWorld's activation, resource and time into its private backup fields for RK2 halfstep.
"""
function rk_backup!( watt::WattWorld)
	watt.backup_a[:] = (p->p.a).(watt.players)
	watt.backup_R    = watt.R
	watt.backup_t    = watt.t
end

#-----------------------------------------------------------------------------------------
"""
	rk_restore!( watt::WattWorld)

Restore WattWorld's activation, resource and time from its private backup fields for RK2 halfstep.
"""
function rk_restore!( watt::WattWorld)
	for (i,p) in pairs(watt.players)
		p.a = watt.backup_a[i]
	end
	watt.R = watt.backup_R
	watt.t = watt.backup_t
end

#-----------------------------------------------------------------------------------------
"""
	omega!( watt::WattWorld, Ω::Float64)

Set WattWorld's omega point (at which stability = 1) equal to Ω.
"""
function omega!( watt::WattWorld, Ω::Float64)
	watt.omega = Ω
end

#-----------------------------------------------------------------------------------------
"""
	feed!( watt::WattWorld, f::Float64)

Set WattWorld's current resource feed value equal to f.
"""
function feed!( watt::WattWorld, f::Float64)
	watt.feed = f
end

#-----------------------------------------------------------------------------------------
"""
	time_generator( watt)

Calculate the time-generator for WattWorld development.
"""
function time_generator( watt::WattWorld)
	apply_regime!(watt)									# Set feed/omega according to current regime
	acts = (p->p.a).(watt.players)						# Activations of Players
	Ks = (p->p.B.K).(watt.players)						# Half-saturation constants K of Players

	activation = (K->saturation(watt.R,K)).(Ks)			# R-dependent activation of Players
	inhibition = saturation(action(watt),-1)			# Players' action inhibits own activation
	net_activation = inhibition * activation			# Net growth of each Player's activation

	da = (net_activation .- BETA_A) .* acts				# Growth - depletion of activations
	dR = watt.feed - BETA_R*watt.R - sum(Ks.*acts)		# Feed rate minus depletion and
														# consumption/production
	(da,dR)
end

#-----------------------------------------------------------------------------------------
"""
	step!( watt::WattWorld, Δt)

Develop WattWorld through a single Runge-Kutta-2 time-step Δt.
"""
function step!( watt::WattWorld, Δt)
	# RK2 half-step:
	rk_backup!(watt)
	Δt2 = Δt/2.0
	da,dR = time_generator( watt)
	for (i,p) in pairs(watt.players)
		p.a += da[i] * Δt2
	end
	watt.R += dR * Δt2
	watt.t += Δt2

	# RK2 full-step:
	da,dR = time_generator( watt)
	rk_restore!(watt)
	for (i,p) in pairs(watt.players)
		p.a += da[i] * Δt
	end
	watt.R += dR * Δt
	watt.t += round(Δt,digits=3)
end

#-----------------------------------------------------------------------------------------
"""
	trajectory( watt::WattWorld, T; Δt=T/500) → Snapshot

Develop WattWorld forward from the current instant through given time duration T.
"""
function trajectory( watt::WattWorld, T::Real; Δt=1.0)
	num_ts = Int(round(T/Δt))+1									# Allow for roundng errors in t
	tfinal = round(watt.t+T,digits=3)
	ts = range(watt.t,tfinal,length=num_ts)
	Δt = step(ts)

	snap = Snapshot( watt)										# Initialise snapshot storage
	snaps = similar([snap],num_ts)
	snaps[1] = snap

	for i in 2:num_ts
		vary!(watt)												# Vary watt world if it is unstable
		step!(watt,Δt)											# Take watt forwards 1 step

		# Repair any invalid values:
		watt.R = max(0.0,min(R_UPB,watt.R))
		for (j,p) in pairs(watt.players)
			if p.a < 0 || p.a > A_UPB
				watt.players[j] = Player()
			end
		end
	
		snaps[i] = Snapshot( watt)								# Store new Snapshot
	end

	return snaps
end

#-----------------------------------------------------------------------------------------
"""
	vary!( watt)

With probability given by WattWorld's current instability, vary the behavioural parameters K and n
of all WattWorld Players by up to an amount DELTA within the respective validity domain.
"""
function vary!( watt::WattWorld)
	instability = 1 - stability(watt)
	deltas = DELTA * instability * (2rand(2,length(watt.players)).-1)
	for (i,p) in pairs(watt.players)
		p.B = (
			K = wrap(p.B.K + deltas[1,i]),
			n = n_LWB + rem( n_SPAN + p.B.n - n_LWB + deltas[2,i], n_SPAN)
		)
	end
end

#-----------------------------------------------------------------------------------------
"""
	stability( watt::WattWorld)

Determine the stability of WattWorld determined by the current resource level R. Currently,
stability is the closeness of R to WattWorld's omega value within a capture well of width
CAPTURE_RADIUS and with monomial degree CAPTURE_CONCAVITY.
"""
function stability( watt::WattWorld)
	1.0 - min( 1.0, (abs(watt.R-watt.omega)/CAPTURE_RADIUS)^CAPTURE_CONCAVITY)
end

#-----------------------------------------------------------------------------------------
"""
	apply_regime!( watt::WattWorld)

Apply WattWorld's current feed/omega test regime for benchmarking.
"""
function apply_regime!( watt::WattWorld)
	if watt.regime <= 1
		# Constant feed and omega:
		watt.omega = 1.0
		watt.feed = 0.2
	elseif watt.regime == 2
		# Single step feed:
		watt.omega = 1.0
		watt.feed = iseven(2watt.t÷DURATION) ? 0.2 : 1.75
	elseif watt.regime == 3
		# Single step omega:
		watt.omega = iseven(2watt.t÷DURATION) ? 1.0 : 2.0
		watt.feed = 0.2
	elseif watt.regime == 4
		# Slow periodic feed:
		watt.omega = 1.0
		watt.feed = iseven(watt.t÷1e6) ? 0.2 : 1.75
	elseif watt.regime == 5
		# Slow periodic omega:
		watt.omega = iseven(watt.t÷1e6) ? 1.0 : 2.0
		watt.feed = 0.2
	else
		# Fast periodic feed with slow periodic omega:
		watt.omega = iseven(watt.t÷8e5) ? 1.0 : 2.0
		watt.feed = iseven(watt.t÷2e5) ? 0.2 : 1.75
	end
end

#-----------------------------------------------------------------------------------------
"""
	regime_string( watt::WattWorld)

Return a short string description of WattWorld's current regime.
"""
function regime_string( watt::WattWorld)
	REGIMES[watt.regime]
end

#-----------------------------------------------------------------------------------------
"""
	report( watt::WattWorld)

Report current status of the WattWorld.
"""
function report( watt::WattWorld)
	println(
		"Time=$(round(watt.t,digits=3)), Resource=$(watt.R), Action=$(action(watt)), ",
		"Stability=$(stability(watt)):"
	)
	for (i,p) in pairs(watt.players)
		if p.a > A_REPORT
			println( "    Player $i: B=$(p.B), a=$(p.a)")
		end
	end
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
"""
	demo()

Build and run WattWorld.
"""
function demo( regime::Int=1)
	watt = WattWorld(N_PLAYERS,regime=regime)
	println(
		"\n======= Simulating a $N_PLAYERS-player WattWorld using regime $regime: ",
		"$(regime_string(watt)) ======="
	)
	report(watt)
	snapshots = trajectory(watt,DURATION,Δt=RK2_STEP)
	report(watt)

	# Display resource, activation and behaviour parameters graphically:
	fig = Figure(fontsize=30,linewidth=5)
	compressed_t = 1:GRAPHICS_COMPRESSION:length(snapshots)
	t_axis = ((s->s.t).(snapshots))[compressed_t]
	ax_R = Axis(fig[1,1], xlabel="time", title="Regime $regime: Global behaviour")
	lines!( ax_R, t_axis, ((s->s.R).(snapshots))[compressed_t], label="Resource")
	lines!( ax_R, t_axis, ((s->s.Ω).(snapshots))[compressed_t], label="Omega")
	lines!( ax_R, t_axis, ((s->s.F).(snapshots))[compressed_t], label="Feed rate")
	Legend( fig[1,2], ax_R)

	ax_a = Axis(fig[1,3], xlabel="time", title="Activation (a)")
	ax_K = Axis(fig[2,1], xlabel="time", title="Half-saturation (B.K)")
	ax_n = Axis(fig[2,3], xlabel="time", title="Cooperativity (B.n)")

	for i in 1:N_PLAYERS
		acts = ((s->s.players[i].a).(snapshots))[compressed_t]
		if any( acts[end÷2:end] .> A_REPORT)
			# Display players with significant activation in the last 1/4 of the simulation:
			lines!( ax_a, t_axis, acts, label="Player $(i)")
			lines!( ax_K, t_axis, ((s->s.players[i].B.K).(snapshots))[compressed_t], 
				label="Player $(i)"
			)
			lines!( ax_n, t_axis, ((s->s.players[i].B.n).(snapshots))[compressed_t],
				label="Player $(i)"
			)
		end
	end
	Legend( fig[2,2], ax_a)
	display( fig)
	fig
end

#-----------------------------------------------------------------------------------------
"""
	publish()

Build and run WattWorld, and save results to the publishable graphics used in
Constructivist Foundations article.
"""
function publish( regime::Int=1)
	watt = WattWorld(N_PLAYERS,regime=regime)
	println(
		"\n======= Simulating a $N_PLAYERS-player WattWorld using regime $regime: ",
		"$(regime_string(watt)) ======="
	)
	report(watt)
	snapshots = trajectory(watt,DURATION,Δt=RK2_STEP)
	report(watt)

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