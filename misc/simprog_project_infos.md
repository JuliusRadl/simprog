# Code

## Project
- Copy Stub.jl as a start

## Structure of ABM (Schelling)
-   Uses Agents Package
-   decleares an agent type Person, derived from the abstract type GridAgent{2}
-   Person agents inherit the fields id::Int and pos::Tuple{2} from GridAgent{2}
-   We define Persons to possess the additional fields comfort::Bool (how comfortable the
    Person feels) and tribe::Int (to which tribe the Person belongs)
-   The method schelling() initialises a Schelling world in which Persons live
-   The method agent_step!() defines how Persons act within their world
-   The method demo() organises, runs and presents the results of the ABM
- Space is grid
```julia
# Setting up a model
ecosys = StandardABM(
        Turtle, # agent struct
        ContinuousSpace(extent); # space
        agent_step!, 
        model_step!,  
        properties # properties dict (StandardABM overwrites getproperty to get these properties
    )
```
- agent_step!() is called once on each agent per step. can be ommited, then 
agent actions must be defined in model_step!()
- model_step!() is called once. 

## Model
- size is in extent as a tuple, eg (60, 60)
```julia
# randomly distribute something in the space (property)
extent = (60, 60)
:algae => rand(Bool, extent)
```
- model properties (like background algae for exampel) are simulated with model_step!()

## Agent
- agents are added with add_agent! (defined by Agents package)
- agents are contained in dict in abm
- continous agents also store velocity (direction in which to move and face ==
a unit vector where its length represents its speed and its angle represents the
direction)
    ```julia
    # 2nd part it anonymous function being passed argument (anonymous function)(argument)
    vel = ecosys.v0 * rand() * (θ->[cos(θ),sin(θ)])(2π*rand())
    # rotate agent by random angle delta up to 12 degrees (pi/15)
    Δ = (2rand()-1)*pi/15
    R = [cos(Δ) -sin(Δ)
            sin(Δ)  cos(Δ)]
    vel = R * collect(vel)
    # More compact:
    cs,sn = (x->(cos(x),sin(x)))((2rand()-1)*pi/15)
    vel = [cs -sn;sn cs]*collect(vel)
    ```
- Properties can be added in struct constructor in a dictionary 
- (there is a distinction between the properties of a struct and the property
:properties of an abm model! abm models overwrite getproperty() to target the
property :properties instead of the real struct properties)
- moved with move_agent!()
- removed with remove_agent!()
- simulated with agent_step!() (contains moving and removing etc.)
```julia
# all they do is in this style:
action!(agent, model)
```
- add new properties to agents
    - define them in the struct
        ```julia
        @agent struct Boid(ContinuousAgent{2,Float64})
            speed::Float64					# Speed of this Boid
        end
        ```
    - add them to the model
        ```julia
        function fields(;
            extent	= (50, 50),		# Dimensions of the world
            dt		= 0.1,			# Time step
            dcAMP	= 1.0,			# Deposition rate of cAMP
            ΛcAMP	= 0.0,			# Diffusion rate of cAMP
            αcAMP	= 0.0,			# Evaporation rate of cAMP
            preference = 0.001,
            aversion = false	
        )
        # Create the Fields space, properties and model:
        model = StandardABM( SlimeMould, ContinuousSpace(extent; spacing=1.0);
            agent_step!, model_step!,
            properties = Dict(		# Model properties applying to all SlimeMoulds:
                :dt			=> dt,
                :dcAMP		=> dcAMP,
                :ΛcAMP		=> ΛcAMP,
                :αcAMP		=> αcAMP,
                :cAMP		=> zeros(Float64, extent),	# Heatmap of cAMP at each location in world
                :preference => preference,
                :aversion => aversion
            )
        )
        ```
    - add them to add_agent!()
        ```julia
        for _ in 1:5
            theta = 2π*rand()								# Random angle
            add_agent!( model, (cos(theta), sin(theta)), preference, aversion)
        end
        ```
    - add them as sliders to your demo
        ```julia
        params = Dict(			
            :ΛcAMP	=> 0:0.01:0.1,
            :αcAMP	=> 0:0.001:0.01,
            :preference => 0:0.001:1.0,
            :aversion => 0:1
        )
        ```


## Initialise Agent
```julia
peter = Schelling.Person(id=1, pos=(1,1), tribe=rand([1, 2]))
fieldtypes(Schelling.Person) # only works with Types, not instantiated objects
```

## General simulation structure
- Initialise ABM
- Fill it with agents
- adapt agent_step!() (gets called for each agent each simulation step)
- adapt run!()

## agent_step!()
- Methods to find nearby agents are implemented in Agents
- Check and update values like comfort etc.
- rules to 
- Move agent with move_agent_single!(agent, model)

## move_agent!()

## run!()
- Data is collect using a vector of tuples where tuple[0] is a symbol (property
of an agent) or a function that will be applied to each agent, and tuple[1] is
a function for aggregating that data into a single value):
```julia
adata = [(:comfort, sum), (agent->agent.pos[1], maximum)]
# Better: Define functions epxlicitly, so their name will be the adata column header
xpos(agent) = agent.pos[1]
adata = [(:comfort, sum), (xpos, maximum)]
# Dangerous: adata is an iterator, a function like this won't work:
partial_mean(data) = mean(data[data .!= 0])
```

## abmvideo
```julia
# import glmakie
using GLMakie
# Define agent colors
tribecolor(agent) = agent.tribe==1 ? :blue : :orange
# Define agent markers
tribemarker(agent) = agent.tribe==1 ? :circle : :rect
# set filename, abm model, marker size, framerate and frames (number of simulation
# timesteps)
abmvideo( "schelling.mp4", schelling(preference);
    title = "Segregation model (comfort threshold $preference)",
    agent_size=10, agent_color=tribecolor, agent_marker=tribemarker,
    framerate = 4, frames = 50
)
```

## Exploratory Playground
```julia
# import palfreymans agent tools
using AgentTools.jl
# Define appearance as in abmvideo above
# first return value is the interactive figure
# second return value is an Observable that can be used to collect data
playground, _ = abmplayground( schelling;
    agent_size = 10, agent_color = tribecolor, agent_marker = tribemarker,
    # set slider for parameters ("Independent Variables") as dict with stepranges
    params = Dict( :preference => 0.0:0.01:1.0)
)
display(playground)
```
- params are for sliders and plotkwargs for plotting
```julia
# Example:
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
```
- get colors from https://juliagraphics.github.io/Colors.jl/stable/colormapsandcolorscales/
or from AgentTools (multicoloured())
- ABMObersable is used to regenerate plots automatically on step (it's the second
return value of abmplayground)

## Field Diffusion, Gradients
- circshift() is used to shift a matrix in a specific direction
```julia
# shifts matrix f by one in the second dimension (columns) (first dimension = rows)
F = circshift(F, (0, 1))
# also works in single dimension (here its a regular column vector)
S = circshift(S, 1)
```
- AgentTools.gradient() applies central difference scheme (numeric approximation
of derivative) to find gradient (get field values in each cardinal direction, then
subtract top from bottom and left from right to get direction in which it increases
the most)
- **Diffraction rate lambda** tells us how much every cell distributes of its field
quantity to each neighboring cell and also receives that fraction from each of its
neighbors

# Theoretical

## Generative Science:
1.  We start from a Research Question about the causes of some collective behaviour;
2.  Our existing theories suggest a Research Hypothesis (HR) about how simple, local agent
    interactions might generate that collective behaviour;
3.  We develop our HR into an Alternative Hypothesis (H1) that describes how we think
    Changes in agent interaction should influence Changes in collective behaviour;
4.  We formulate a Null Hypothesis (H0) that Denies the influence proposed by H1;
5.  We construct and perform an experiment to Disprove H0;
6.  If the experiment successfully disproves H0, this justifies our belief in H1 and HR.

## GS applied by Schelling:
1.  His research question: To what extent does cultural discomfort drive segregation?
2.  HR: Agents generate segregation by relocating to reduce cultural discomfort;
3.  H1: If HR is correct (that is, if discomfort-reduction causes segregation), then
    Decreasing agents' level of comfort should Increase their level of segregation.
4.  H0: Decreasing agents' comfort level makes segregation decrease or stay constant.
5.  Our model must allow agents to relocate on the basis of their preference for tribally
    similar neighbours, and we must be able to vary this individual preference and
    measure the corresponding level of segregation in the community.

## Escape Plateaus
- Schelling Agents get stuck because they can't move freely once minimum viable comfort
is achieved
- In evolution plateaus are escaped by mutation == randomness

## Non-Computability & Embodiment
- difficult to predict the community's behaviour based on the actions of its individual agents
- Idee von Embodiment: Intelligenz / Denken setzt (physische) Verkörperung voraus =>
ist nicht nur Software, die auch auf Computer laufen kann
- Embodied Systems haben einen berechenbaren und einen unberechenbaren / chaotischen Teil

## Supervenience / Complexity / Emergence
- the behavioural identity of a living organism **supervenes** on the physical structure of its body; that is, we cannot compute the identity from the body, but if two bodies are identical, then their identities must also be identical
- Identity **emerges** from the organism's physical structure
- **Complexity** is the ability of a system to generate collective, **non-computable** behaviour that is **Reliable**, but different each time
- if collective, non-computable behaviour is also reliable, we describe it as being **Complex**.

## Decisions / Choosing
- Decisions are not made by individuals or environement, but by their interaction (ecosystem)


## Reference Mode
- Reference Mode: Einfacheres Verhalten, gesehen in realen Systemen, dass das ABM-System zeigen können
muss, bevor die Forschungsfrage beantwortet werden kann / um als gültiges System durchzugehen
- Was ist “sustainably periodic collective behavior? Wichtig ist, dass das Verhalten nachhaltig
gezeigt wird. **Mit den richtigen Einstellungen** schwingt die Population auch => periodisch
- Welches kollektive Verhalten der Schildkröten können wir hier beobachten, welches
aus den Interaktionen individueller Agenten entsteht? Dass die Zahl der Agenten und Algen konstant bleibt
- Was muss gelten, damit Embodiment gewährleistet ist? Das Verhalten muss kollektiv sein, und nicht vorgegeben

## Fields
- are their own entity
- provide useful information to organisms that they use to guide their behavior
- interaction between agents and fields creates niche

# Stabilisation
- Provided some non-local field (such as the twig) is present that can link together the
influence of small, randomly occurring fluctuations in the environment, this field can
consolidate these fluctuations into larger scale stability.