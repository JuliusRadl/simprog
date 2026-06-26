#========================================================================================#
"""
    Ecosystem

To-do: Specify the Ecosystem module.

Author: Niall Palfreyman (January 2025)
"""
module Ecosystem
include( "../../Development/Generative/AgentTools.jl")
using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
    Turtle

To-do: Define the Turtle Agent type.
"""
@agent struct Turtle(ContinuousAgent{2,Float64})
    energy::Float64				# My current energy
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
    ecosystem( max_speed=5; extent=(60,60))

To-do: Initialise the Ecosystem model.
"""
function ecosystem(;
    dt          = 0.1,          # Time-step interval for the model
    v0          = 5.0,          # Maximum initial speed of Turtles
    E0          = 20,           # Maximum initial energy of Turtles
    extent      = (60,60)       # Spatial extent of the model
)
    # To-do: Initialise the model properties
    properties = Dict(
        :dt     => dt,
        :v0     => v0,
        :E0     => E0,
    )
    ecosys = StandardABM(
        Turtle,
        ContinuousSpace(extent);
        agent_step!,
        properties
    )

    # To-do: Initialise the agents
    n_agents = 1
    for _ in 1:n_agents
        vel = ecosys.v0 * [1,1]
        energy = ecosys.E0
        add_agent!( ecosys; vel, energy)
    end

    ecosys
end

#-----------------------------------------------------------------------------------------
"""
    agent_step!( Turtle, model)

To-do: Define the Turtles' behaviour.
"""
function agent_step!( me::Turtle, model)
    # To-do: Set state

    # To-do: Perceive

    # To-do: Decide

    # To-do: Act

    # To-do: Move
    move_agent!( me, model, model.dt)

    return
end

#-----------------------------------------------------------------------------------------
"""
    demo()

Demonstrate the Ecosystem model.
"""
function demo()
    params = Dict(
        # To-do: Specify model exploration parameters
        :v0             => 0:0.1:10,
        :E0             => 1:50,
    )
    plotkwargs = (
        # To-do: Specify plotting keyword arguments
        agent_size      = 10,
        agent_color     = :red,
        agent_marker    = :circle,
        adata           = [(:energy, sum)],
        mdata           = [nagents],
    )
    playground, _ = abmplayground( ecosystem; params, plotkwargs...)
    display(playground)
end

end # ... of module Ecosystem