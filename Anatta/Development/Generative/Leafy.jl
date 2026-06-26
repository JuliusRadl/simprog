#========================================================================================#
"""
    Leafy (Non-Computability copy)

This module generates collective behaviour structures that display Non-Computability.
Non-computable properties are not necessarily irreducible, since they are certainly determined by
the individual behaviours of system components. Nevertheless, non-computable properties are
Reliable yet non-predictable. That is, we can generate them reliably out of randomised individual
beviours, yet we cannot predict them using any computational process based purely on those
individual behaviours.

Author: Niall Palfreyman, February 2025.
"""
module Leafy
include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
    Builder

A Builder is a very simple Agent with no internal state other than the default properties id,
pos, vel. The only thing we really need from a Builder in this simulation is its position.
"""
@agent struct Builder(ContinuousAgent{2,Float64})
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
    noncomputability()

Initialise the NonComputability model.
"""
function noncomputability(;
    n_builders	= 1,            # Number of Builders building the non-computable structure
    r			    = 0.5,          # Length of Builders' jump towards next food source
    n_sources	    = 3,            # Number of sources of of food in the model
    extent		    = (15,15),      # Spatial extent of the model
    base_point      = [7.5, 7.5]        # Initial Point
)
    properties = Dict(
        :n_builders => n_builders,
        :r 				=> r,
        :sources		=> [],
        :n_sources		=> n_sources,
        :spu			=> 30,
        :footprints		=> Vector{Point2f}(undef,0),
        :base_point     => base_point,
        :affine         => [
            (P -> ([0 0; 0 0.16] * P + [0,0])),
            (P -> ([0.85 0.04; -0.04 0.85] * P + [0,1.6])),
            (P -> ([0.20 -0.26; 0.23 0.22] * P + [0,1.6])),
            (P -> ([-0.15 0.28; 0.26 0.24] * P + [0,0.44]))
        ]
    )

    model = StandardABM( Builder,
        ContinuousSpace(extent);
        agent_step!, model_step!,
        properties
    )

    # Create regular polygonal distribution of food sources:
    Δθ = 2π/n_sources
    worldwidth = minimum(extent)
    model.sources = map(0:Δθ:2π-Δθ) do θ
        Tuple(0.5*worldwidth*([-sin(θ),cos(θ)] .+ 1))
    end

    # To-do: Initialise the agents
    for _ in 1:n_builders
        # Create Builders at random initial positions in the world:
        add_agent!( Tuple(base_point), model, vel=[0.0,0.0])
    end
    
    model
end

#-----------------------------------------------------------------------------------------
"""
    agent_step!( builder, model)

Builder jumps a constant fraction r of the distance towards a randomly selected source.
"""
function agent_step!( builder, model)
    p = collect(builder.pos) - model.base_point
    # generate random number between 0-1
    r = rand()
    # check if number is in intervall ]0-0.01]
    if (r < 0.01)
        p = model.affine[1](p)
    elseif (r < 0.08)
        p = model.affine[3](p)
    elseif (r < 0.15)
        p = model.affine[4](p)
    else
        p = model.affine[2](p)
    end
    move_agent!( builder, Tuple(p+model.base_point), model)
end

"""
    model_step!(model)

After all agents have moved one step, record their current position as a footprint.
"""
function model_step!(model)
    append!( model.footprints, [Point2f(p.pos) for p in allagents(model)])
end

#-----------------------------------------------------------------------------------------
"""
    sources(model)

Return all source locations of the NonComputability model.
"""
function sources(model)
    Point2f.(model.sources)
end

#-----------------------------------------------------------------------------------------
"""
    footprints(model)

Return all footprints of the NonComputability model.
"""
function footprints(model)
    model.footprints
end

#-----------------------------------------------------------------------------------------
"""
    demo()

Demonstrate the NonComputability model.
"""
function demo()
    params = Dict(
        # To do: Playground slider values:
        :n_builders	=> 1:10,
        :r				=> 0:0.1:1,
        :n_sources		=> 3:7,
        :spu		    => 30:100
    )
    plotkwargs = (
        # To-do: Specify plotting keyword arguments
        agent_size      = 20,
        agent_color     = :green,
        agent_marker    = :circle,
    )

    # Create playground displaying all Builders and sources:
    playground, abmobs = abmplayground( noncomputability; params, plotkwargs...)
    # Add (Observable) sources and footprints
    scatter!( lift( sources, abmobs.model), color=:blue,  markersize=20)
    scatter!( lift( footprints, abmobs.model), color=:red, markersize=1)

    display(playground)
end

end # ... of module NonComputability
