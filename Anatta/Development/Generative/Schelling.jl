#========================================================================================#
"""
    Schelling

This module reproduces Schelling's model of segregation in social communities.

Author: Niall Palfreyman (December 2024)
"""
module Schelling

# To-do: Add graphics backend package
using Agents

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
    Person

A Person is a very simple Agent existing in a discrete 2-dimensional grid space. The @agent macro
ensures that the Person type autonmatically inherits the implicit state variables `id` and
`position`, and we explicitly add here the state variables `comfort` and `tribe`.
"""
@agent struct Person(GridAgent{2})
    comfort::Bool = false		            # All Persons start life feeling uncomfortable
    tribe::Int					            # Each Person's tribe (0/1) will need initialising
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
    schelling( preference=1.0; extent=(60,60))

Initialise the size, stepping function, properties and initial population of a Schelling model.
"""
function schelling( preference=1.0; extent=(60,60))
    # To-do: Define preference as a model property
    properties = Dict(
    )

    schelling_model = StandardABM(
        Person,
        GridSpace(extent);
        agent_step!,
        properties
    )

    # To-do: Set n_agents to 80% of the total number of grid points
    n_agents = 1
    for n in 1:n_agents
        # To-do: Place Persons of random tribe
        add_agent_single!( schelling_model; tribe=0)
    end

    return schelling_model
end

#-----------------------------------------------------------------------------------------
"""
    agent_step!( Person, model)

Define the evolution rule acting once per step on each activated Person agent:
    If my neighbourhood group contains at least a certain proportion of members of my own tribe,
    I feel comfortable; otherwise, I feel uncomfortable and react accordingly.
"""
function agent_step!( me::Person, model)
    # Calculate the proportion of neighbours belonging to my tribe:
    n_nbrs = n_mytribe = 0
    for nbr in nearby_agents(me,model)
        n_nbrs += 1
        if me.tribe == nbr.tribe
            n_mytribe += 1
        end
    end
    proportion_mytribe = (n_nbrs > 0) ? n_mytribe/n_nbrs : 0.0

    # To-do: Decide how to react
    me.comfort = proportion_mytribe ≥ 1.0

    # To-do: If uncomfortable, jump to a random empty grid location

    return
end

#-----------------------------------------------------------------------------------------
"""
    present_insight()

Calculate the average results of several trials of the Schelling model for each value of
`preference` in the range 0.0:0.2:1.0. For efficiency, halt each trial when all agents are happy.
Present these results as a graph demonstrating that rising individual preference for tribally
similar neighbours drives community segregation. This graph plots preference along the horizontal
axis and segregation along the vertical axis. We define the level of segregation of a community as
the proportion of agents that possess a maximal number (8 in the Chebyshev metric) of tribally
similar neighbours.
"""
function present_insight()
    preferences = 0.0:0.02:1.0
    segregation = zeros(length(preferences))
    ntrials = 5
    halt_condition(abm,t) = all(map(agent->agent.comfort,allagents(abm))) || t>99

    for (i,preference) in enumerate(preferences)
        for _ in 1:ntrials
            abm = schelling(preference)
            step!( abm, halt_condition)
            for agent in allagents(abm)
                nbrs = nearby_agents(agent,abm)
                n_tribal_nbrs = 0
                # To-do: Count tribally similar neighbours

                if n_tribal_nbrs == 8
                    segregation[i] += 1             # This agent is fully surrounded by own tribe
                end
            end
        end
    end
    segregation = segregation ./ (ntrials*nagents(schelling()))
    lines( preferences, segregation, axis = (
        title="Individual preference drives community segregation",
        xlabel="Preference", ylabel="Segregation"
    ))
end

#-----------------------------------------------------------------------------------------
"""
    demo()

Demonstrate the Schelling model.
"""
function demo(preference=1.0)
    # Cap preference values above and below:
    preference = max(0.0,min(1.0,preference))
    abm = schelling(preference)

    # To-do: Define tribe position data
    adata = [(:comfort, sum)]
    agent_df, model_df = run!( abm, 9; adata)
    agent_df

    # To-do: Generate video output

    # To-do: Create an exploratory playground

end

end # ... of module Schelling