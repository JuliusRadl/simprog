#========================================================================================#
"""
    Ecosystem

To-do: Specify the Ecosystem module.

Author: Niall Palfreyman (January 2025)
"""
module Ecosystem

include("../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
    Turtle

To-do: Define the Turtle Agent type.
"""
@agent struct Turtle(ContinuousAgent{2,Float64})
    energy::Float64# My current energy
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
    ecosystem( max_speed=5; extent=(60,60))

To-do: Initialise the Ecosystem model.
"""
function ecosystem(;
    dt=0.1,             # Time-step interval for the model
    v0=5.0,             # Maximum initial speed of Turtles
    E0=20,              # Maximum initial energy of Turtles
    extent=(60, 60)     # Spatial extent of the model
)
    # Done: Initialise the model properties
    properties = Dict(
        :dt => 0.1,     # Time-step interval for the model
        :n_turtles => 5,       # Initial number of turtles
        :v0 => 5,       # Maximum initial speed of a turtle
        :E0 => 100.0,   # Maximum initial energy of a turtle
        :Δliving => 1.0,     # Energy cost of living
        :Δeating => 7.0,      # Energy benefit of eating one algae
        :algae => rand(Bool,extent), # Random algae grid in the background
        :prob_regrowth  => 0.01    # Prob. that any algae will regrow in a step
    )
    ecosys = StandardABM(
        Turtle,
        ContinuousSpace(extent);
        agent_step!,
        model_step!,
        properties
    )

    # Done: Initialise the agents
    n_agents = 5
    for _ in 1:n_agents
        vel = ecosys.v0 * rand() * (θ->[cos(θ),sin(θ)])(2π*rand())
        energy = rand(1:ecosys.E0)
        add_agent!(ecosys; vel, energy)
    end

    ecosys
end

#-----------------------------------------------------------------------------------------
"""
    agent_step!( Turtle, model)

To-do: Define the Turtles' behaviour.
"""
function agent_step!(me::Turtle, model)
    # To-do: Set state
    me.energy -= model.Δliving * model.dt
    if me.energy < 0
        remove_agent!( me, model)
        return
    end

    # To-do: Perceive

    # To-do: Decide

    # To-do: Act
    eat!( me, model)
    reproduce!( me, model)

    # To-do: Move
    cs,sn = (x->(cos(x),sin(x)))((2rand()-1)*pi/15)
    me.vel = [cs -sn;sn cs]*collect(me.vel)
    move_agent!( me, model, model.dt)

    return
end

#-----------------------------------------------------------------------------------------
"""
    eat!( turtle, model)

Causes a turtle agent to eat / remove an algae from its current position in
the model.
"""

function eat!( turtle, model)
        indices = get_spatial_index( turtle.pos, model.algae, model)
        if model.algae[indices]
            turtle.energy += model.Δeating
            model.algae[indices] = false
        end
    end

#-----------------------------------------------------------------------------------------
"""
    reproduce!( parent::Turtle, model)

Creates a new Turtle if 
"""
function reproduce!( parent::Turtle, model)
        if parent.energy > model.E0 && rand() < 0.01
            parent.energy -= model.E0
            add_agent!( parent.pos, model, parent.vel, model.E0)
    end
end

#-----------------------------------------------------------------------------------------
"""
    model_step!( model)

Advance the model, causing algae to grow etc.
"""
function model_step!( model)
        empty_locs = .!model.algae
        model.algae[empty_locs] .= (rand(count(empty_locs)).<model.prob_regrowth)
end

#-----------------------------------------------------------------------------------------
"""
    demo()

Demonstrate the Ecosystem model.
"""
function demo()
    params = Dict(
        :prob_regrowth  => 0:0.0001:0.01,
        :E0                 => 10.0:200.0,
        :Δeating        => 0:0.1:10.0,
        :Δliving        => 0:0.1:10.0,
    )
    plotkwargs = (
        # To-do: Specify plotting keyword arguments
        agent_size=10,
        agent_color=multicoloured,
        agent_marker=wedge,
        adata=[(a->isa(a,Turtle),count)], alabels=["Turtles"],
        mdata=[(m->sum(m.algae))], mlabels=["Algae"],
        heatarray       = (model->model.algae),
        heatkwargs      = (colormap=[:black,:darkgreen],colorrange=(0,1)),
        add_colorbar    = false
    )
    playground, _ = abmplayground(ecosystem; params, plotkwargs...)
    display(playground)
end

end # ... of module Ecosystem