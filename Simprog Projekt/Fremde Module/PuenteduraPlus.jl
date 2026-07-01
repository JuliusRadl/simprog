#=================================================================================================#
"""
	PuenteduraPlus

Eine Simulation, welche die genetische und entwicklungsbedingte Suche implementiert.
Diese ist aufbauend auf Puenteduras Ergebnisse und wird erweitert, um zu erforschen, wie sich die 
verschiedene Entwicklungsschritte T auf die Suchdauer verhalten.

Author: Charlotte Weinmüller 28.06.2025

Forschungsfrage: Welche mathematische Funktion beschreibt die Dauer (bis 10% Agenten das Zielgenom 
                gefunden haben) im Mittel für verschiedene, steigende T-Werte?

Forschungshypothese: In Abhängigkeit von steigenden T-Werten sinkt die mittlere Dauer streng 
                    monoton auf einen minimalen Grenzwert > 0, da Agenten T-Anzahl an 
                    Enwicklungsschritten pro Simulationsschritt machen können.

Nullhypothese: Die Kurve der mittleren Dauer ist nicht monoton sinkend.
	
"""
module PuenteduraPlus

include("./AgentTools.jl")

using Agents, GLMakie, .AgentTools 

#--------------------------------------------------------------------------------------------------
# Modultyp:
#--------------------------------------------------------------------------------------------------
"""
	DevAgent

DevAgents haben eigene Genome und versuchen durch Reproduktion mit Crossover und Mutation und
durch Fitness das Zielgenom zu finden. Abhängig von ihrer Fitness resultiert die Energy und Farbe.
"""
@agent struct DevAgent(ContinuousAgent{2,Float64})
    genome::Vector{Bool}        # Festes Genom
    p_genome::Vector{Bool}      # Genom für Plastizität (plastisch wenn Allel = 1)
    fitness::Float64            # Fitness des Agenten
    energy::Float64             # Lebensenergie des Agenten, basierend auf der Fitness
    colour::Symbol              # Farbgebung, um die Fitness der Agenten zu visualisieren
end

#--------------------------------------------------------------------------------------------------
# Modul Methoden:
#--------------------------------------------------------------------------------------------------
"""
	puentedura_plus()

Initialisierung der PuenteduraPlus Simulation.
"""
function puentedura_plus(; 
    crossover=true,
    mu_rate=1e-3,
    genome_length=20,
    T=1000,
    development=true,
    elitism = -1
)
    extent = (40,40)
    target = fill(true, genome_length)  # Zielgenom

    properties = Dict(
        :crossover      => crossover,           # Erfolgt bei der Reproduktion ein Crossover?
        :mu_rate        => mu_rate,             # Wie hoch ist die Wahrscheinlichkeit zu mutieren?
        :target         => target,              # Zielgenom
        :T              => T,                   # Entwicklungsschritte
        :mean           => 1.0,                 # mittlere Fitness der Population 
		:sigma			=> 0.0,			        # Standardabweichung der Fitness
        :max_pop        => 0.2*prod(extent),    # Maximal erlaubte Populationsgröße
        :development    => development,         # Sind Agenten entwicklungsfähig?
        :life_energy    => 10.0,                # Lebensenergie für den Start ins Leben
        :living_cost    => 0.1,                 # Wie viel Lebensenergie kostet das Leben?
        :elitism        => elitism,             # Elitismus
        :bestID         => 1,                   # ID des Agenten mit der besten Fitness 
    )

    # ABM-Modell definieren
    model = StandardABM(DevAgent, ContinuousSpace(extent);
        model_step!, agent_step!, properties)

    # Agenten mit dem Positionen und Genomen zuweisen
    for i in 1:model.max_pop/2  
        pos = (rand()*extent[1], rand()*extent[2])
        genome = rand([false, true], genome_length)
        p_genome = rand([false, true], genome_length)
        
        # Agenten ins ABM hinzufügen
        baby_dev = add_agent!(model, pos, genome, p_genome, 1.0, 0.0, :white)  

        # Jeweilige Energie, die zum reproduzieren und überleben gebraucht wird, berechnen
        set_energy!(baby_dev, model)
    end

    return model
end


#--------------------------------------------------------------------------------------------------
"""
	agent_step!( dev, model)

Ein Schritt im Leben eines DevAgenten des PuenteduraPlus Versuchs.
Der Agent bewegt sich und verliert dadurch Energie.
- Bei genügend Energie: Reproduktion (Crossover + Mutation)
- Bei zu wenig Energie: Agent stirbt und wird aus dem model entfernt
"""
function agent_step!(dev::DevAgent, model)
    walk!(dev, model)
    if dev.energy < 0
        remove_agent!(dev, model)
        return
    end
    reproduce(dev, model)
end


#--------------------------------------------------------------------------------------------------
"""
	model_step!( model)

Globale Berechnungen: für den aktuellen Mittelwert (mean) und Standardabweichung.
Und Ermitteln des Agenten mit der höchsten Fitness aus der gesamten Population.
"""
function model_step!(model)
    alldevs = collect(allagents(model))

    # prüfen ob die Population noch existiert
    if isempty(alldevs)
        println("Population ausgestorben...")
        return
    end

    fitnesses = [dev.fitness for dev in alldevs]
    model.mean = mean(fitnesses)
    model.sigma = std(fitnesses)
    model.bestID = alldevs[findmax(fitnesses)[2]].id
end

#--------------------------------------------------------------------------------------------------
"""
	set_energy!( dev, model)

Berechnet die Energie eines Agenten basierend auf individueller und gesamter Populations-Fitness.
Ist die eigene Fitness überdurchschnittlich, steigt die Energie; bei Energie < 0 stirbt der Agent.
Die Energy beeinflusst Überleben und Reproduktion.
"""
function set_energy!(dev::DevAgent, model)
    # Fitness für den jeweiligen Agenten berechnen
    dev.fitness = set_fitness!(dev.genome, model.target, model, dev.p_genome)

    if model.sigma==0
		# Wenn alle devs gleich fit sind: 
		energy = 1.0
	else
        energy = (dev.fitness - model.mean)/(model.sigma) - model.elitism 
    end
    dev.energy = energy < 0 ? -1 : (energy * model.life_energy)

    # Farben des Agenten, bedingt durch seine Fitness, setzen
    if dev.fitness == 1.0
		dev.colour = :green
	elseif dev.fitness > 0.9
		dev.colour = :yellow
	elseif dev.fitness > 0.5
		dev.colour = :red
	elseif dev.fitness > 0.2
		dev.colour = :gray
	else
		dev.colour = :black
	end

end


#--------------------------------------------------------------------------------------------------
"""
	set_fitness!( genome, target)

Berechnet die Fitness eines Genoms im Vergleich zum Target.
Für einen modulareren Austausch an Fitnessfunktionen werden ausgelagerte Funktionen aufgerufen.
"""
function set_fitness!( genome::Vector{Bool}, target::Vector{Bool}, model, p_genome::Vector{Bool})
    # Ohne Entwicklung
    if model.development == false 
        fit = fit_without_dev(genome, target)
    end

    # Mit Entwicklung 
    if model.development == true
        fit = fit_with_dev(genome, target, model, p_genome)
    end

    return fit
end

#--------------------------------------------------------------------------------------------------
"""
	fit_without_dev(genome, target)

Fitnessberechnung ohne Entwicklung: Fitness = 1.0 bei exaktem Treffer, sonst 0.0
"""
function fit_without_dev(genome::Vector{Bool}, target::Vector{Bool})
    fit = genome == target ? 1.0 : 0.0

    return fit
end

#--------------------------------------------------------------------------------------------------
"""
	fit_with_dev(genome, target, model, p_genome)

Fitnessberechnung mit Entwicklung: Bis zu T Versuche durch plastische Änderungen.
Die erste erfolgreiche Runde bestimmt die Hinton-Nowlan-Fitness, normiert auf [0, 1].
"""
function fit_with_dev(genome::Vector{Bool}, target::Vector{Bool}, model, p_genome::Vector{Bool})
    T = model.T
    len = length(genome)

    # Wenn kein Entwicklungsschritt: aus Effizienzgründen Wert direkt ohne Berechnung zurückgeben
    if T == 0
        return fit = genome == target ? 1.0 : 0.0
    end

    phenotype = copy(genome)                # kopieren, da Phenotypen nicht genetisch fixiert sind
    found_round = -1                        # zum initialisieren auf unerreichbaren Wert setzen

    # Entwicklungsversuche (n = 0 bis T)
    for n in 0:T
        # Prüfen, ob phenotype schon exakt dem target entspricht
        if phenotype == target
            found_round = n
            break
        end
        # Plastische Stellen zufällig verändern
        rand_vec = rand(Bool, len)          # 50% der Stellen sollen mit rand(Bool) bekommen 
        to_change = rand_vec .& p_genome    # wenn beides true, dann ist Entwicklung möglich
        phenotype[to_change] .= rand(Bool)  # Wo to_change == true: neuer boolescher Wert
    end

    # Falls in keiner der T Runden getroffen, dann found_round = T (worst-case)
    if found_round == -1
        found_round = T
    end

    # Hinton-Nowlan-Fitness: je geringer found_round ist, desto höher die Fitness
    hn_fit = 1 + ((len - 1) * (T - found_round)) / T
    fit = (hn_fit - 1.0) / (len - 1)        # normiert für mehr Flexibilität
    return fit    
end

#--------------------------------------------------------------------------------------------------
"""
	walk!( dev, model)

Der Agent bewegt sich und verliert währenddessen Energy zum Leben.
"""
function walk!( dev::DevAgent, model)
	wiggle!(dev)
	move_agent!(dev, model, rand())
	dev.energy -= model.living_cost * model.life_energy
end


#--------------------------------------------------------------------------------------------------
"""
	reproduce(dev, model)

Eltern bekommen bei ausreichend Fitness zwei Kinder mit ähnlichen Attributen, die aber
etwas verändert sein können (durch Mutation und Crossover).
"""
function reproduce(mummy::DevAgent, model)
    len = length(mummy.genome)
    # Energy, die die Reproduktion kostet
    birth_cost = model.living_cost * model.life_energy
    # Reproduktion, wenn das model noch nicht maximal voll ist und die Energy ausreicht
    if nagents(model) < model.max_pop && mummy.energy > birth_cost
        daddy = random_nearby_agent(mummy, model)
        if daddy !== nothing && daddy.energy > birth_cost
            # Start mit Crossover
            billy_genome, billy_p_genome, sally_genome, sally_p_genome = 
            crossover(mummy.genome, daddy.genome, mummy.p_genome, daddy.p_genome, model)

            # ... dann mutieren ...
            billy_genome, billy_p_genome, sally_genome, sally_p_genome = 
            mutate!(len, billy_genome, billy_p_genome, sally_genome, sally_p_genome, model)

            # ... und die beiden Kinder in die Welt setzen ...
            pos1 = (rand() * 40, rand() * 40)
            pos2 = (rand() * 40, rand() * 40)
            set_energy!(add_agent!(
                model, pos1, billy_genome, billy_p_genome, 1.0, 0.0, :white), model)
            set_energy!(add_agent!(
                model, pos2, sally_genome, sally_p_genome, 1.0, 0.0, :white), model)

            # ... und Energy der Eltern um Geburtskosten verringern
            mummy.energy -= birth_cost
            daddy.energy -= birth_cost
        end
    end
end

#--------------------------------------------------------------------------------------------------
"""
	crossover(genom_mum, genom_dad, p_genome_mum, p_genome_dad, model)

Bei der Reproduktion werden Genome der Eltern an je einem zufälligen Crossoverpunkt rekombiniert, 
getrennt für genome und p_genome. So kann ein Allel unabhängig fix oder lernfähig sein.
"""
function crossover(genom_mum::Vector{Bool}, genom_dad::Vector{Bool}, 
        p_genome_mum::Vector{Bool}, p_genome_dad::Vector{Bool}, model)

    len = length(genom_mum)
    # Crossover Punkte (fixe Gen und das plastische Gen unabhängig voneinander random wählen)
    xpt_fixed = model.crossover ? rand(1:len-1) : len
    xpt_plastic = model.crossover ? rand(1:len-1) : len

    # Erzeuge für Kind1 (Billy) die Genome fix & plastisch
    billy_genome = vcat(genom_dad[1:xpt_fixed], genom_mum[xpt_fixed+1:end])
    billy_p_genome = vcat(p_genome_dad[1:xpt_plastic], p_genome_mum[xpt_plastic+1:end])	
    
    # Erzeuge für Kind2 (Sally) die jeweiligen Genome 
    sally_genome = vcat(genom_mum[1:xpt_fixed], genom_dad[xpt_fixed+1:end])
    sally_p_genome = vcat(p_genome_mum[1:xpt_plastic], p_genome_dad[xpt_plastic+1:end])

    # Neue Genome zurückgeben
    return billy_genome, billy_p_genome, sally_genome, sally_p_genome
end


#--------------------------------------------------------------------------------------------------
"""
	mutate!(len, billy_genome, billy_p_genome, sally_genome, sally_p_genome, model)
In die Funktion wird die Länge des Genoms, die vier Genome der beiden Kinder und das ABM-Model übergeben.
Dann zufällige Mutation der Genome der Kinder.
"""
function mutate!(len::Int, billy_genome, billy_p_genome, sally_genome, sally_p_genome, model)
    # Billys Genom - fix & plastisch - mutieren
    mutated_loci = rand(len) .< model.mu_rate
    billy_genome[mutated_loci] = rand([true,false],count(mutated_loci))
    mutated_loci = rand(len) .< model.mu_rate
    billy_p_genome[mutated_loci] = rand([true,false],count(mutated_loci))

    # Sallys Genom - fix & plastisch - mutieren
    mutated_loci = rand(len) .< model.mu_rate
    sally_genome[mutated_loci] = rand([true,false],count(mutated_loci))
    mutated_loci = rand(len) .< model.mu_rate
    sally_p_genome[mutated_loci] = rand([true,false],count(mutated_loci))

    # Mutierte Genome zurückgeben
    return billy_genome, billy_p_genome, sally_genome, sally_p_genome
end


#--------------------------------------------------------------------------------------------------
"""
	calculate_mean_duration(T)

Simuliert für jeden T-Wert repeats (100) Mal die Evolution bis 10 % der Agenten Fitness 1.0 
erreichen oder die maximale Schrittzahl überschritten ist. 
Gibt die durchschnittliche Dauer zurück.
"""
function calculate_mean_duration(T; repeats::Int=100, maximum_steps = 100000)
    # Vektor zur Speicherung der Schritte pro Wiederholung initialiseren
    times = Vector{Float64}() 

        for i in 1:repeats
            println("Nun startet die $i-te Wiederholung")
            # Initialisiere das Modell mit dem aktuellen T-Wert
            model = puentedura_plus(T=T)
            # Führe die Simulation durch, bis die Fitness von 1.0 oder maximum_steps erreicht
            step_count = 0
            while step_count < maximum_steps
                step!(model, 1) # Führe einen Schritt der Simulation durch
                step_count += 1
                # Überprüfe, ob die Fitness von 1.0 erreicht wurde
                all_agents = allagents(model)
                successful_agents = count(a -> a.fitness == 1.0, all_agents)
                total_agents = length(all_agents)
                # Abbruchbedingung: mindestens 10% erfolgreiche Agenten
                if total_agents > 0 && successful_agents >= 0.1 * total_agents
                    println("Ich war erfolgreich im $step_count")
                    break
                end
                if isempty(allagents(model))
                    println("Die Population ist bei $T im Schritt $step_count ausgestorben")
                    break
                end
            end
            push!(times, step_count)
        end

    # mittlere Dauer für den T-Wert zurückgeben
    return mean(times)
end


#--------------------------------------------------------------------------------------------------
"""
	plot_graph!()

Simuliert die mittlere Evolutionsdauer für verschiedene Werte von T (Entwicklungsdauer) und 
visualisiert den Einfluss auf die Such-Effizienz. Die Resultate werden geplottet und gespeichert.
"""
function plot_graph!(; T_values=0:50:2000, repeats::Int=100, maximum_steps::Int=100000)
    T_vals = collect(T_values)          # Liste der T-Werte
    n = length(T_vals)
    means = Vector{Float64}(undef, n)   # Speicher für die Mittelwerte der Dauer
    
    # Simuliere mittlere Suchdauer für jeden T-Wert
    for (idx, T) in enumerate(T_vals)
        println("Starte Simulationen für T = $T ...")
        mean_duration = calculate_mean_duration(T;repeats=repeats, maximum_steps=maximum_steps)
        means[idx] = mean_duration
    end

    # Plot initialisieren
    fig = Figure(resolution = (800, 600))   # Erstelle ein neues Figure-Objekt
    ax = Axis(fig[1, 1],                    # Füge eine Achse hinzu
        xlabel = "T",
        ylabel = "Mean Duration",
        title = "Influence of T on Search Efficiency"
    )
    # Linie zeichnen
    lines!(ax, T_vals, means)
    # Punkte markieren
    scatter!(ax, T_vals, means)

    display(fig)

    # Plot als Datei speichern
    save("plot_dev.png", fig)
end



#-----------------------------------------------------------------------------------------
"""
	demo()

Demonstriere den PuenteduraPlus Versuch.
"""
function demo()
    model = puentedura_plus()

    # Parameter für die interaktive Steuerung
    params = Dict(
		:crossover  => false:true,
		:mu_rate    => 0:0.001:0.1,
		:development => false:true,
        :T          => 0:5:2000,
        :elitism    => -5:1:0
	)
	plotkwargs = (
		agent_color=(dev->dev.colour),
		agent_size=20,
		mdata = [(m->m.mean), (m->m[m.bestID].fitness)],
		mlabels = ["Mean Fitness", "Best Fitness"],
	)

    # Starte interaktive Simulation mit Parametern und Plotoptionen
	playground,abmplt = abmplayground( puentedura_plus; params, plotkwargs...)

    # Dynamisches Beobachten des besten Genoms und seiner Plastizität
    best_genome = lift((wld -> join(Int.(wld[wld.bestID].genome))), abmplt.model)
    best_p_genome = lift((wld -> join(Int.(wld[wld.bestID].p_genome))), abmplt.model)

    # Als Text in Plot anzeigen
	text!( 10,1.1, text=best_genome, color=:red, fontsize=16, align=(:left,:top))
    text!( 10,1.05, text=best_p_genome, color=:blue, fontsize=16, align=(:left,:top))

	playground
end

end     # ... of module PuenteduraPlus 