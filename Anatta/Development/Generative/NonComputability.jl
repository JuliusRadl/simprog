#========================================================================================#
"""
    Non-Computability

This module generates collective behaviour structures that display Non-Computability.
Non-computable properties are not necessarily irreducible, since they are certainly determined by
the individual behaviours of system components. Nevertheless, non-computable properties are
Reliable yet non-predictable. That is, we can generate them reliably out of randomised individual
beviours, yet we cannot predict them using any computational process based purely on those
individual behaviours.

Author: Niall Palfreyman, February 2025.
"""
module NonComputability
include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
    Grasshopper

A Grasshopper is a very simple Agent with no internal state other than the default properties id,
pos, vel. The only thing we really need from a Grasshopper in this simulation is its position.
"""
@agent struct Grasshopper(ContinuousAgent{2,Float64})
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
    noncomputability()

Initialise the NonComputability model.
"""
function noncomputability(;
    n_grasshoppers	= 1,            # Number of Grasshoppers building the non-computable structure
    r			    = 0.5,          # Length of Grasshoppers' jump towards next food source
    n_sources	    = 4,            # Number of sources of of food in the model
    extent		    = (100,100)     # Spatial extent of the model
)
    properties = Dict(
        :n_grasshoppers => n_grasshoppers,
        :r 				=> r,
        :sources		=> [],
        :n_sources		=> n_sources,
        :spu			=> 30,
        :footprints		=> Vector{Point2f}(undef,0)
    )

    noncompute = StandardABM( Grasshopper,
        ContinuousSpace(extent);
        agent_step!, model_step!,
        properties
    )

    # Create regular polygonal distribution of food sources:
    Δθ = 2π/n_sources
    worldwidth = minimum(extent)
    noncompute.sources = map(0:Δθ:2π-Δθ) do θ
        Tuple(0.5*worldwidth*([-sin(θ),cos(θ)] .+ 1))
    end

    # To-do: Initialise the agents
    for _ in 1:n_grasshoppers
        # Create grasshoppers at random initial positions in the world:
        add_agent!( Tuple(worldwidth*rand(2)), noncompute, vel=[0.0,0.0])
    end
    
    noncompute
end

#-----------------------------------------------------------------------------------------
"""
    agent_step!( grasshopper, noncompute)

Grasshopper jumps a constant fraction r of the distance towards a randomly selected source.
"""
function agent_step!( grasshopper, noncompute)
    source = rand(noncompute.sources)                                   # Select a random source
    move_agent!( grasshopper,                                           # Jump towards selected src
        ((1-noncompute.r).*grasshopper.pos .+ noncompute.r.*source),    # using convex combination
        noncompute
    )
end

"""
    model_step!(noncompute)

After all agents have moved one step, record their current position as a footprint.
"""
function model_step!(noncompute)
    append!( noncompute.footprints, [Point2f(p.pos) for p in allagents(noncompute)])
end

#-----------------------------------------------------------------------------------------
"""
    sources(noncompute)

Return all source locations of the NonComputability model.
"""
function sources(noncompute)
    Point2f.(noncompute.sources)
end

#-----------------------------------------------------------------------------------------
"""
    footprints(noncompute)

Return all footprints of the NonComputability model.
"""
function footprints(noncompute)
    noncompute.footprints
end

#-----------------------------------------------------------------------------------------
"""
    demo()

Demonstrate the NonComputability model.
"""
function demo()
    params = Dict(
        # To do: Playground slider values:
        :n_grasshoppers	=> 1:10,
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

    # Create playground displaying all grasshoppers and sources:
    playground, abmobs = abmplayground( noncomputability; params, plotkwargs...)
    # Add (Observable) sources and footprints
    scatter!( lift( sources, abmobs.model), color=:blue,  markersize=20)
    scatter!( lift( footprints, abmobs.model), color=:red, markersize=1)

    display(playground)
end

end # ... of module NonComputability
