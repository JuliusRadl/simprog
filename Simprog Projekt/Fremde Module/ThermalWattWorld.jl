#=========================================================================================#
# spatial thermal regulation through narrative stabilisation
#
# The model separates fast activation, slower structural variation, and movement so that
# their respective contributions to local thermal stability can be measured.
#=========================================================================================#
"""
    ThermalWattWorld

A spatial, temperature-based extension of WattWorld.

Agents regulate a diffusing temperature field through short-term activations and mutable
heating or cooling parameters. Agents interact only through this shared field: local thermal
success reduces both structural variation and movement speed. Test regimes introduce cold,
heat, or repeated pulses.
"""
module ThermalWattWorld

using Agents
using GLMakie
using Random
using Statistics

const agenttools_candidates = [
    joinpath(@__DIR__, "AgentTools.jl"),
    joinpath(@__DIR__, "AgentTools(1).jl"),
]

const agenttools_path = let index = findfirst(isfile, agenttools_candidates)
    isnothing(index) && error("Place AgentTools.jl in the same directory as ThermalWattWorld.jl.")
    agenttools_candidates[index]
end

include(agenttools_path)
using .AgentTools

#=========================================================================================#
# agents and model construction
#
# Each agent combines a structural regulator K with a fast activation state. Negative K
# values heat a local cell; positive K values cool it.
#=========================================================================================#
"""
    ThermalAgent

A mobile local regulator. `K` is a mutable heating or cooling parameter, while
`activation` is the short-term strength with which the agent currently acts.
"""
@agent struct ThermalAgent(ContinuousAgent{2, Float64})
    K::Float64
    activation::Float64
end

"""
    disturbance_events(test_regime, cold_amplitude, heat_amplitude)

Return a reproducible schedule of ambient-temperature disturbances. `:baseline` contains
no disturbance and is the default experimental regime.
"""
function disturbance_events(test_regime, cold_amplitude, heat_amplitude)
    events = NamedTuple{(:start, :duration, :offset), Tuple{Int, Int, Float64}}[]

    if test_regime == :baseline
        return events
    elseif test_regime == :cold_pulse
        push!(events, (start = 1000, duration = 100, offset = -cold_amplitude))
    elseif test_regime == :heat_pulse
        push!(events, (start = 1000, duration = 100, offset = heat_amplitude))
    elseif test_regime == :repeated_cold
        push!(events, (start = 500, duration = 100, offset = -cold_amplitude))
        push!(events, (start = 1000, duration = 100, offset = -cold_amplitude))
        push!(events, (start = 1500, duration = 100, offset = -cold_amplitude))
    elseif test_regime == :alternating
        push!(events, (start = 500, duration = 100, offset = -cold_amplitude))
        push!(events, (start = 1000, duration = 100, offset = heat_amplitude))
        push!(events, (start = 1500, duration = 100, offset = -cold_amplitude))
    else
        error(
            "Unknown test_regime $(repr(test_regime)). Use :baseline, :cold_pulse, " *
            ":heat_pulse, :repeated_cold, or :alternating.",
        )
    end

    return events
end

"""
    thermal_watt_world(; kwargs...)

Create a two-dimensional temperature field populated by mobile heating and cooling agents.

Set `stabilise = true` for local K stabilisation. `stabilise_movement` defaults to the
same value, but can be set separately to isolate structural and behavioural mechanisms.
"""
function thermal_watt_world(;
    extent = (50, 50),
    n_agents = 200,
    dt = 0.1,
    target_temperature = 15.0,
    ambient_temperature = 10.0,
    diffusion_rate = 0.25,
    relaxation_rate = 0.04,
    thermal_strength = 3.00,
    activation_decay = 0.02,
    activation_maximum = 3.0,
    mutation_rate = 0.06,
    K_minimum = 0.10,
    K_maximum = 2.00,
    variation_tolerance = 3.0,
    movement_tolerance = 10.0,
    minimum_speed = 0.001,
    maximum_speed = 0.5,
    turning_angle = pi,
    cold_amplitude = 8.0,
    heat_amplitude = 8.0,
    test_regime = :baseline,
    stabilise = true,
    stabilise_movement = stabilise,
    seed = nothing,
)
    !isnothing(seed) && Random.seed!(seed)

    temperature = fill(ambient_temperature, extent)
    thermal_input = zeros(Float64, extent)
    events = disturbance_events(test_regime, cold_amplitude, heat_amplitude)

    properties = Dict(
        :dt => dt,
        :temperature => temperature,
        :thermal_input => thermal_input,
        :target_temperature => target_temperature,
        :ambient_temperature => ambient_temperature,
        :diffusion_rate => diffusion_rate,
        :relaxation_rate => relaxation_rate,
        :thermal_strength => thermal_strength,
        :activation_decay => activation_decay,
        :activation_maximum => activation_maximum,
        :mutation_rate => mutation_rate,
        :K_minimum => K_minimum,
        :K_maximum => K_maximum,
        :variation_tolerance => variation_tolerance,
        :movement_tolerance => movement_tolerance,
        :minimum_speed => minimum_speed,
        :maximum_speed => maximum_speed,
        :turning_angle => turning_angle,
        :test_regime => test_regime,
        :stabilise_K => stabilise,
        :stabilise_movement => stabilise_movement,
        :events => events,
        :step_count => 0,
        :last_ambient_temperature => ambient_temperature,
        :mean_K_change => 0.0,
        :mean_movement_speed => 0.0,
        :temperature_volatility => 0.0,
        :K_change_total => 0.0,
        :speed_total => 0.0,
    )

    model = StandardABM(
        ThermalAgent,
        ContinuousSpace(extent; spacing = 1.0);
        agent_step!,
        model_step!,
        properties,
    )

    for _ in 1:n_agents
        angle = 2.0*pi*rand()
        velocity = (cos(angle), sin(angle))
        add_agent!(model, velocity, random_K(model), 1.00)
    end

    return model
end

"""
    random_K(model)

Draw one initial structural regulator. The gap around zero avoids agents whose thermal
effect is too weak to make their regulatory role observable.
"""
function random_K(model)
    magnitude = model.K_minimum + rand()*(model.K_maximum - model.K_minimum)

    return rand(Bool) ? magnitude : -magnitude
end

#=========================================================================================#
# local WattWorld-like regulation
#
# The local temperature substitutes for WattWorld's common resource. Activation responds
# quickly to each agent's local temperature, while K becomes sticky under successful
# regulation. Diffusion alone provides coupling between agents.
#=========================================================================================#
"""
    hill_response(signal, K)

Return a signed Hill-style response. Cooling agents respond more strongly to larger
signals, while heating agents respond more strongly to smaller signals.
"""
function hill_response(signal, K)
    magnitude = abs(K)
    rising_response = signal/(signal + magnitude)

    return K > 0.0 ? rising_response : 1.0 - rising_response
end

"""
    local_temperature(agent, model)

Return the temperature at the grid cell currently occupied by `agent`.
"""
function local_temperature(agent, model)
    index = get_spatial_index(agent.pos, model.temperature, model)

    return model.temperature[index]
end

"""
    activation_change(agent, model, temperature)

Calculate the fast activation change from the local temperature and the agent's own K.
Agents remain coupled indirectly through the shared diffusing temperature field.
"""
function activation_change(agent, model, temperature)
    signal = max(0.0, temperature)/model.target_temperature
    response = hill_response(signal, agent.K)

    return agent.activation*(response - model.activation_decay)
end

"""
    variation_scale(agent, model, temperature)

Return the permitted K variation. Under local stabilisation, well-regulated agents become
structurally sticky; the control condition preserves a constant variation rate.
"""
function variation_scale(agent, model, temperature)
    !model.stabilise_K && return model.mutation_rate

    local_error = abs(temperature - model.target_temperature)

    return model.mutation_rate*min(1.0, local_error/model.variation_tolerance)
end

"""
    mutate_K!(agent, model, temperature)

Apply bounded stochastic structural variation. Crossing the gap around zero changes an
agent's regulatory role only when variation is sufficiently large to pass the gap.
"""
function mutate_K!(agent, model, temperature)
    previous_K = agent.K
    candidate = previous_K + variation_scale(agent, model, temperature)*randn()
    candidate = clamp(candidate, -model.K_maximum, model.K_maximum)

    if abs(candidate) < model.K_minimum
        direction = candidate == 0.0 ? sign(previous_K) : sign(candidate)
        candidate = direction*model.K_minimum
    end

    agent.K = candidate

    return abs(agent.K - previous_K)
end

#=========================================================================================#
# movement and field dynamics
#
# Agents remain mobile, but movement is slower in locally successful thermal conditions.
# This makes behavioural exploration and structural variation comparable stabilisation paths.
#=========================================================================================#
"""
    movement_speed(model, temperature)

Return a random-walk speed that increases with local temperature error. Movement
stabilisation can be disabled independently to isolate the effect of structural K change.
"""
function movement_speed(model, temperature)
    !model.stabilise_movement && return model.maximum_speed

    local_error = abs(temperature - model.target_temperature)
    error_fraction = min(1.0, local_error/model.movement_tolerance)

    return model.minimum_speed +
        error_fraction*(model.maximum_speed - model.minimum_speed)
end

"""
    move_randomly!(agent, model, temperature)

Turn by a random angle and move according to the temperature-dependent exploration speed.
"""
function move_randomly!(agent, model, temperature)
    speed = movement_speed(model, temperature)

    wiggle!(agent, model.turning_angle)
    move_agent!(agent, model, speed)

    return speed
end

"""
    agent_step!(agent, model)

Update fast activation, add local heating or cooling, vary K, and perform a random walk.
Thermal effects are accumulated before the field update to avoid dependence on agent order.
"""
function agent_step!(agent, model)
    temperature = local_temperature(agent, model)
    activation_delta = activation_change(agent, model, temperature)
    agent.activation = clamp(
        agent.activation + model.dt*activation_delta,
        0.0,
        model.activation_maximum,
    )

    index = get_spatial_index(agent.pos, model.temperature, model)
    model.thermal_input[index] += -model.thermal_strength*agent.activation*agent.K
    model.K_change_total += mutate_K!(agent, model, temperature)
    model.speed_total += move_randomly!(agent, model, temperature)
end

"""
    ambient_at_step(model, step)

Return the ambient temperature specified by the active test regime at one simulation step.
"""
function ambient_at_step(model, step)
    offset = 0.0

    for event in model.events
        active = event.start <= step < event.start + event.duration
        active && (offset += event.offset)
    end

    return model.ambient_temperature + offset
end

"""
    model_step!(model)

Diffuse and relax the field, apply accumulated agent effects, and finalise the per-step
metrics used to evaluate thermal performance and dynamical stability.
"""
function model_step!(model)
    next_step = model.step_count + 1
    previous_temperature = copy(model.temperature)
    ambient = ambient_at_step(model, next_step)
    diffused = diffuse4(model.temperature, model.diffusion_rate*model.dt)
    relaxed = diffused .+ model.dt*model.relaxation_rate .* (ambient .- diffused)
    updated = relaxed .+ model.dt .* model.thermal_input

    model.temperature[:] = clamp.(updated, 0.0, 40.0)
    model.temperature_volatility = Statistics.mean(
        abs.(model.temperature .- previous_temperature),
    )
    model.mean_K_change = model.K_change_total/max(1, nagents(model))
    model.mean_movement_speed = model.speed_total/max(1, nagents(model))
    model.last_ambient_temperature = ambient
    model.step_count = next_step

    fill!(model.thermal_input, 0.0)
    model.K_change_total = 0.0
    model.speed_total = 0.0
end

#=========================================================================================#
# performance metrics and test routines
#
# Field error measures accuracy, habitable coverage measures spatial success, and
# volatility measures whether success persists rather than appearing only momentarily.
#=========================================================================================#
"""
    mean_temperature_error(model)

Return the mean absolute error of the full temperature field relative to the target.
"""
function mean_temperature_error(model)
    return Statistics.mean(abs.(model.temperature .- model.target_temperature))
end

"""
    mean_agent_temperature_error(model)

Return the mean local temperature error experienced at the positions of all agents.
"""
function mean_agent_temperature_error(model)
    errors = (
        abs(local_temperature(agent, model) - model.target_temperature)
        for agent in allagents(model)
    )

    return Statistics.mean(errors)
end

"""
    habitable_fraction(model; band = 2.0)

Return the fraction of field cells within `band` degrees of the target temperature.
"""
function habitable_fraction(model; band = 2.0)
    within_band = abs.(model.temperature .- model.target_temperature) .<= band

    return count(within_band)/length(within_band)
end

"""
    stability_index(model; band = 2.0)

Combine habitable coverage with temporal field stability. Report its components alongside
this index because a single score can conceal different failure modes.
"""
function stability_index(model; band = 2.0)
    coverage = habitable_fraction(model; band)
    normalised_volatility = model.temperature_volatility/band

    return coverage/(1.0 + normalised_volatility)
end

"""
    temperature_volatility(model)

Return the mean absolute temperature change of the full field during the latest step.
"""
function temperature_volatility(model)
    return model.temperature_volatility
end

"""
    performance_metrics(model)

Return the current observable measures for one simulation step.
"""
function performance_metrics(model)
    return (
        step = model.step_count,
        ambient_temperature = model.last_ambient_temperature,
        field_error = mean_temperature_error(model),
        agent_error = mean_agent_temperature_error(model),
        habitable_fraction = habitable_fraction(model),
        stability_index = stability_index(model),
        temperature_volatility = model.temperature_volatility,
        mean_K_change = model.mean_K_change,
        mean_movement_speed = model.mean_movement_speed,
    )
end

"""
    recovery_time(records, event; baseline_window = 80, margin = 0.10)

Return the number of steps after one disturbance until field error returns near its
pre-disturbance average. Return `missing` when recovery is not observed in the run.
"""
function recovery_time(records, event; baseline_window = 80, margin = 0.10)
    baseline_start = max(1, event.start - baseline_window)
    baseline_errors = [
        record.field_error
        for record in records
        if baseline_start <= record.step < event.start
    ]

    isempty(baseline_errors) && return missing

    threshold = Statistics.mean(baseline_errors)*(1.0 + margin)
    recovery_start = event.start + event.duration
    recovered_index = findfirst(
        record -> record.step >= recovery_start && record.field_error <= threshold,
        records,
    )

    isnothing(recovered_index) && return missing

    return records[recovered_index].step - recovery_start
end

"""
    summarise_test(records, events)

Return final performance and recovery measures for one complete run.
"""
function summarise_test(records, events)
    final_window_start = max(1, length(records) - 99)
    final_records = records[final_window_start:end]
    recovery_times = [recovery_time(records, event) for event in events]
    observed_recoveries = [time for time in recovery_times if !ismissing(time)]

    return (
        final_field_error = Statistics.mean(record.field_error for record in final_records),
        final_agent_error = Statistics.mean(record.agent_error for record in final_records),
        final_habitable_fraction = Statistics.mean(
            record.habitable_fraction for record in final_records
        ),
        final_stability_index = Statistics.mean(
            record.stability_index for record in final_records
        ),
        final_temperature_volatility = Statistics.mean(
            record.temperature_volatility for record in final_records
        ),
        final_movement_speed = Statistics.mean(
            record.mean_movement_speed for record in final_records
        ),
        recovery_times = recovery_times,
        mean_recovery_time = isempty(observed_recoveries) ? missing :
            Statistics.mean(observed_recoveries),
    )
end

"""
    run_test(; steps = 1400, stabilise = true, seed = 1, kwargs...)

Run one condition and return time-series records, the disturbance schedule, and a compact
summary. Use a non-baseline `test_regime` to evaluate recovery after temperature shocks.
"""
function run_test(; steps = 1400, stabilise = true, seed = 1, kwargs...)
    model = thermal_watt_world(; stabilise, seed, kwargs...)
    records = NamedTuple[]

    for _ in 1:steps
        step!(model, 1)
        push!(records, performance_metrics(model))
    end

    return (
        records = records,
        events = model.events,
        summary = summarise_test(records, model.events),
        stabilise = stabilise,
    )
end

"""
    compare_conditions(; n_runs = 20, steps = 1400, test_regime = :repeated_cold, kwargs...)

Run matched stabilisation and control trials over multiple seeds. This supplies the raw
per-run results needed for a later statistical comparison rather than relying on one run.
"""
function compare_conditions(;
    n_runs = 20,
    steps = 1400,
    test_regime = :repeated_cold,
    kwargs...,
)
    results = NamedTuple[]

    for seed in 1:n_runs
        stabilised = run_test(
            ;
            steps,
            stabilise = true,
            seed,
            test_regime,
            kwargs...,
        )
        control = run_test(
            ;
            steps,
            stabilise = false,
            seed,
            test_regime,
            kwargs...,
        )

        for result in (stabilised, control)
            summary = result.summary
            push!(results, (
                seed = seed,
                condition = result.stabilise ? "stabilisation" : "control",
                final_field_error = summary.final_field_error,
                final_agent_error = summary.final_agent_error,
                final_habitable_fraction = summary.final_habitable_fraction,
                final_stability_index = summary.final_stability_index,
                final_temperature_volatility = summary.final_temperature_volatility,
                final_movement_speed = summary.final_movement_speed,
                mean_recovery_time = summary.mean_recovery_time,
            ))
        end
    end

    return results
end

#=========================================================================================#
# interactive exploration
#
# The default playground has no disturbance. Pass a named test regime to inspect how the
# same model behaves under cold, heat, repeated-cold, or alternating perturbations.
#=========================================================================================#
"""
    demo(; test_regime = :baseline, stabilise = true,
         stabilise_movement = stabilise, kwargs...)

Open an Anatta playground.

Supported regimes are `:baseline`, `:cold_pulse`, `:heat_pulse`, `:repeated_cold`, and
`:alternating`. The default is `:baseline`, so no environmental disturbance is applied.
Set `stabilise_movement = false` to test K stabilisation without movement stabilisation.
"""
function demo(;
    test_regime = :baseline,
    stabilise = true,
    stabilise_movement = stabilise,
    kwargs...,
)
    constructor = () -> thermal_watt_world(
        ;
        test_regime,
        stabilise,
        stabilise_movement,
        kwargs...,
    )

    plotkwargs = (
        agent_color = agent -> agent.K < 0.0 ? :firebrick : :deepskyblue,
        agent_marker = :circle,
        agent_size = 10,
        heatarray = model -> model.temperature .- model.target_temperature,
        heatkwargs = (colormap = :bwr, colorrange = (-10.0, 10.0)),
        add_colorbar = true,
        colorbar_label = "Temperature difference T - Ω",
        mdata = [
            mean_temperature_error,
            mean_agent_temperature_error,
            habitable_fraction,
            stability_index,
            temperature_volatility,
            model -> model.step_count,
            model -> model.mean_movement_speed,
        ],
        mlabels = [
            "Field error",
            "Agent error",
            "Habitable fraction",
            "Stability index",
            "Temperature volatility",
            "Simulation step",
            "Mean movement speed",
        ],
    )

    playground, = abmplayground(constructor; plotkwargs...)

    display(playground)

    return playground
end

"""
    smoke_test(; kwargs...)

Execute ten model steps without the graphical interface and return the resulting metrics.
This isolates model execution from Anatta's playground controls when diagnosing a run issue.
"""
function smoke_test(; kwargs...)
    model = thermal_watt_world(; kwargs...)

    step!(model, 10)

    return performance_metrics(model)
end

export ThermalAgent
export compare_conditions
export demo
export habitable_fraction
export mean_agent_temperature_error
export mean_temperature_error
export performance_metrics
export run_test
export stability_index
export smoke_test
export thermal_watt_world
export temperature_volatility

end
