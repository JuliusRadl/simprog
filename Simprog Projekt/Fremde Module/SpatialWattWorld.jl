#========================================================================================#
# spatial temperature regulation through local narrative stabilisation
#========================================================================================#
"""
    SpatialWattWorld

A spatial extension of WattWorld in which agents regulate a diffusing temperature field.
Each agent has a mutable regulator parameter `K` and a rapidly changing activation.
Successful local regulation reduces mutation of `K`; the control condition keeps mutation
independent of local temperature error.
"""
module SpatialWattWorld

using Agents, GLMakie, Random, Statistics

const agenttools_candidates = [
    joinpath(@__DIR__, "AgentTools.jl"),
    joinpath(@__DIR__, "AgentTools(1).jl"),
]

const agenttools_path = let path = findfirst(isfile, agenttools_candidates)
    isnothing(path) && error("Place AgentTools.jl in the same directory as SpatialWattWorld.jl.")
    agenttools_candidates[path]
end

include(agenttools_path)
using .AgentTools

#========================================================================================#
# agent and model definitions
#========================================================================================#
"""
    ThermalAgent

An agent with a mutable heating or cooling parameter `K` and a fast activation state.
Negative `K` heats its local cell; positive `K` cools it.
"""
@agent struct ThermalAgent(ContinuousAgent{2,Float64})
    K::Float64
    activation::Float64
end

"""
    thermal_world(; kwargs...)

Create a two-dimensional temperature field populated by local heating and cooling agents.
Set `stabilise = true` for the narrative-stabilisation condition and `false` for the
control condition.
"""
function thermal_world(;
    extent = (80, 80),
    n_agents = 500,
    dt = 0.1,
    target_temperature = 2.0,
    ambient_temperature = 3.0,
    diffusion_rate = 0.8,
    relaxation_rate = 0.03,
    heating_strength = 1.20,
    activation_decay = 0.25,
    mutation_rate = 0.06,
    tolerance = 3.0,
    disturbance_amplitude = 8.0,
    disturbance_start = 5000,
    disturbance_duration = 100,
    stabilise = true,
    seed = nothing,
)
    !isnothing(seed) && Random.seed!(seed)

    temperature = fill(ambient_temperature, extent)
    properties = Dict(
        :dt => dt,
        :temperature => temperature,
        :target_temperature => target_temperature,
        :ambient_temperature => ambient_temperature,
        :diffusion_rate => diffusion_rate,
        :relaxation_rate => relaxation_rate,
        :heating_strength => heating_strength,
        :activation_decay => activation_decay,
        :mutation_rate => mutation_rate,
        :tolerance => tolerance,
        :disturbance_amplitude => disturbance_amplitude,
        :disturbance_start => disturbance_start,
        :disturbance_duration => disturbance_duration,
        :stabilise => stabilise,
        :step_count => 0,
        :mean_k_change => 0.0,
        :k_change_total => 0.0,
    )

    model = StandardABM(
        ThermalAgent,
        ContinuousSpace(extent; spacing = 1.0);
        agent_step!,
        model_step!,
        properties,
    )

    for _ in 1:n_agents
        K = rand(Bool) ? rand(0.10:0.01:2.00) : -rand(0.10:0.01:2.00)
        add_agent!(model, (cos(2.0*pi*rand()), sin(2.0*pi*rand())), K, 0.10)
    end

    return model
end

#========================================================================================#
# local regulation and structural variation
#========================================================================================#
"""
    hill_response(temperature, K)

Return the activation tendency of a regulator at a local temperature. Heating regulators
with negative `K` respond more strongly to cold cells; cooling regulators with positive
`K` respond more strongly to warm cells.
"""
function hill_response(temperature, K, target_temperature)
    magnitude = abs(K)
    signal = max(0.0, temperature)/target_temperature
    response = signal/(signal + magnitude)

    return K > 0.0 ? response : 1.0 - response
end

"""
    local_temperature(agent, model)

Read the temperature field at the spatial cell occupied by `agent`.
"""
function local_temperature(agent, model)
    index = get_spatial_index(agent.pos, model.temperature, model)

    return model.temperature[index]
end

"""
    activation_change(agent, model, temperature)

Calculate the rapid activation dynamics. The multiplicative term prevents activation from
becoming negative and keeps the response proportional to the current active capacity.
"""
function activation_change(agent, model, temperature)
    response = hill_response(temperature, agent.K, model.target_temperature)
    mean_drive = Statistics.mean(
        abs(other.K)*other.activation for other in allagents(model)
    )
    coordination = 0.8 + mean_drive/(1.0 + mean_drive)

    return agent.activation*(coordination*response - model.activation_decay)
end

"""
    mutation_scale(agent, model, temperature)

Return the permitted magnitude of structural variation. In the stabilisation condition,
local success makes `K` sticky. The control condition removes this success-dependent rule.
"""
function mutation_scale(agent, model, temperature)
    !model.stabilise && return model.mutation_rate

    error = abs(temperature - model.target_temperature)

    return model.mutation_rate*min(1.0, error/model.tolerance)
end

"""
    mutate_K!(agent, model, temperature)

Apply bounded stochastic variation to a regulator parameter. The gap around zero avoids
regulators that contribute almost no heating or cooling effect.
"""
function mutate_K!(agent, model, temperature)
    old_K = agent.K
    scale = mutation_scale(agent, model, temperature)
    candidate = clamp(old_K + scale*randn(), -2.0, 2.0)

    if abs(candidate) < 0.10
        candidate = sign(candidate == 0.0 ? old_K : candidate)*0.10
    end

    agent.K = candidate

    return abs(agent.K - old_K)
end

"""
    agent_step!(agent, model)

Update rapid activation, locally heat or cool the temperature field, and then vary the
slower regulator structure.
"""
function agent_step!(agent, model)
    temperature = local_temperature(agent, model)
    change = activation_change(agent, model, temperature)
    agent.activation = clamp(agent.activation + model.dt*change, 0.0, 3.0)

    index = get_spatial_index(agent.pos, model.temperature, model)
    model.temperature[index] -= model.dt*model.heating_strength*agent.activation*agent.K

    model.k_change_total += mutate_K!(agent, model, temperature)
end

#========================================================================================#
# spatial field dynamics and disturbances
#========================================================================================#
"""
    current_ambient_temperature(model)

Return the environmental temperature at the current step. A finite cold pulse supplies a
repeatable disturbance from which recovery can be measured.
"""
function current_ambient_temperature(model)
    within_pulse = model.disturbance_start <= model.step_count <
        model.disturbance_start + model.disturbance_duration

    return within_pulse ? model.ambient_temperature - model.disturbance_amplitude :
        model.ambient_temperature
end

"""
    model_step!(model)

Diffuse temperature, relax the field toward the current ambient temperature, and finalise
per-step structural-change data after all agents have acted.
"""
function model_step!(model)
    model.step_count += 1
    ambient = current_ambient_temperature(model)
    diffused = diffuse4(model.temperature, model.diffusion_rate*model.dt)
    relaxed = diffused .+ model.dt*model.relaxation_rate .* (ambient .- diffused)

    model.temperature[:] = relaxed[:]
    model.mean_k_change = model.k_change_total/max(1, nagents(model))
    model.k_change_total = 0.0
end

#========================================================================================#
# measurement and interactive exploration
#========================================================================================#
"""
    mean_temperature_error(model)

Return the mean absolute difference between local temperature and the target temperature.
"""
function mean_temperature_error(model)
    return Statistics.mean(abs.(model.temperature .- model.target_temperature))
end

"""
    habitable_fraction(model; band = 2.0)

Return the fraction of cells whose temperature lies within `band` degrees of the target.
"""
function habitable_fraction(model; band = 2.0)
    return Statistics.mean(abs.(model.temperature .- model.target_temperature) .<= band)
end

"""
    run_experiment(; steps = 1200, stabilise = true, seed = 1, kwargs...)

Run one condition and return a table of the temperature error, habitable fraction, and
mean structural variation at every model step.
"""
function run_experiment(; steps = 1200, stabilise = true, seed = 1, kwargs...)
    model = thermal_world(; stabilise, seed, kwargs...)
    records = NamedTuple[]

    for step in 1:steps
        step!(model, 1)
        push!(records, (
            step = step,
            mean_error = mean_temperature_error(model),
            habitable_fraction = habitable_fraction(model),
            mean_k_change = model.mean_k_change,
            condition = stabilise ? "stabilisation" : "control",
        ))
    end

    return records
end

"""
    demo(; stabilise = true)

Open an Anatta playground for one condition. Run this once with `stabilise = true` and
once with `stabilise = false` to inspect the experimental and control dynamics.
"""
function demo(; stabilise = true)
    constructor = () -> thermal_world(; stabilise)
    plotkwargs = (
        agent_color = agent -> agent.K < 0.0 ? :firebrick : :deepskyblue,
        agent_marker = :circle,
        agent_size = 12,
        heatarray = model -> model.temperature,
        heatkwargs = (colormap = :thermal, colorrange = (0.0, 10.0)),
        add_colorbar = true,
        mdata = [mean_temperature_error, habitable_fraction],
        mlabels = ["Mean temperature error", "Habitable fraction"],
    )

    return abmplayground(constructor; plotkwargs...)[1]
end

export ThermalAgent, demo, habitable_fraction, mean_temperature_error, run_experiment,
    thermal_world

end
