#=================================================================================================#
"""
	NarrativeExploration

Investigate the characteristics of Puentedura's implementation of narrative exploration, in which
agents seek solutions not only through selection of genotypes, but also by narratively exploring
different phenotypic expressions of their genomes.

Author: Niall Palfreyman, February, 2026
"""
module NarrativeExploration

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#--------------------------------------------------------------------------------------------------
# Module types:
#--------------------------------------------------------------------------------------------------
"""
	NarrativeAgent

NarrativeAgents possess two genomes: a value genome and a plasticity genome of equal length. The
phenotypic expression of loci at which the plasticity allele is false, is equal to the allele at
the corresponding value genome locus. At loci where the plasticity locus is true, the phenotypic
expression is chosen randomly from the available alleles on each life-step.
"""
@agent struct NarrativeAgent(ContinuousAgent{2,Float64})
	v_genome::Vector			# Value genome
	p_genome::Vector{Bool}		# Plasticity genome (allele true indicates plasticity at this locus)
	phenotype::Vector			# Phenotype resulting from plastic expression of genome
	energy::Float64				# Current energy of the agent
	dissonance::Float64			# Dissonance of agent with its environment
	plasticity::Float64			# Amount of plasticity in agent's genome
	min_diss::Float64			# Best dissonance achieved so far
	colour::Symbol				# Colour indicating agent's energy status
	fixed::Bool					# Fix an agent's phenotype by ceasing all further exploration
end

#--------------------------------------------------------------------------------------------------
# Module constants:
#--------------------------------------------------------------------------------------------------
const benchtest_alphabet = 'A':'Z'

#--------------------------------------------------------------------------------------------------
# Module methods:
#--------------------------------------------------------------------------------------------------
"""
	narrative_exploration()

Initialise the NarrativeExploration Simulation.
"""
function narrative_exploration(;
	narrative = false,
	mu_rate = 2e-3,
	benefit = 2.0,
	difficulty = 1.0,
	allele_count = 2,
	target_length = 21,
	life_energy = 2.0,
	unimodal_fitness = true,
)
	pPop = 0.3											# Probability of a cell becoming populated
	extent = (40,40)
	living_cost = 0.1									# Energy cost of each life-step
	alleles, target = create_problem(allele_count,target_length)
	allele_count = length(alleles)
	target_length = length(target)

	properties = Dict(
		:narrative				=> narrative,			# Using narrative exploration?
		:mu_rate				=> mu_rate,				# Allele mutation probabiliy
		:benefit				=> benefit,				# Benefit of correctly guessing the target
		:unimodal_fitness		=> unimodal_fitness,	# Use unimodal fitness landscape?
		:difficulty				=> difficulty,			# Difficulty of comparing larger dissonances
		:life_energy			=> life_energy,			# Average initial energy of agents
		:living_cost			=> living_cost,			# Energy cost of each life-step
		:birth_cost				=> 2living_cost,		# Energy cost of reproducing
		:search_speed			=> 1.0,					# Agents' speed
		:interaction_distance	=> 3.0,					# Max interaction distance between agents
		:alleles				=> alleles,				# Alphabet of available allele symbols
		:allele_count			=> allele_count,		# Number of available alleles
		:target					=> target,				# Target string to be discovered
		:target_length			=> target_length,		# Length of target string
		:problem_size			=> target_length*(allele_count-1),
		:max_pop				=> pPop*prod(extent),	# Maximum (approx) allowed population size
		:mean_diss				=> 0.5,					# Mean dissonance of population
		:mean_plasticity		=> 0.5,					# Mean plasticity of population
		:min_dissID				=> 1,					# ID of agent with minimum dissonance
	)
	ne = StandardABM( NarrativeAgent, ContinuousSpace(extent);
		model_step!, agent_step!, properties
	)
	genome_length = length(ne.target)
	for i in 1:ne.max_pop/2
		θ = 2π * rand()
		v_genome = rand(ne.alleles,genome_length)
		p_genome = rand(Bool,genome_length)
		natalia = add_agent!( ne, (cos(θ),sin(θ)), v_genome, p_genome, copy(v_genome),
			rand() * ne.life_energy,					# Initial energy
			1.0, 1.0, 1.0,								# Provisional dissonances
			:black,										# Provisional colour
			false										# Not fixed
		)
		set_dissonance!(natalia, ne)
	end

	return ne
end

#--------------------------------------------------------------------------------------------------
"""
	model_step!( ne)

Calculate and save current global mean and minimum dissonance values for entire model population.
"""
function model_step!(ne)
	allnatalias = collect(allagents(ne))

	if isempty(allnatalias)
		println("Population extinct ...")
		return
	end

	ne.mean_diss = mean([natalia.dissonance for natalia in allnatalias])
	ne.mean_plasticity = mean([natalia.plasticity for natalia in allnatalias])
	ne.min_dissID = allnatalias[findmin([natalia.min_diss for natalia in allnatalias])[2]].id
end

#--------------------------------------------------------------------------------------------------
"""
	agent_step!( natalia, ne)

One life-step of a NarrrativeAgent: age, explore development, reproduce and possibly die.
"""
function agent_step!(natalia::NarrativeAgent, ne)
	walk!( natalia, ne)
	explore!( natalia, ne)
	if natalia.energy <= 0
		remove_agent!(natalia, ne)
		return
	end
	reproduce!( natalia, ne)
end

#--------------------------------------------------------------------------------------------------
"""
	walk!( natalia, ne)

Agent moves one step and loses a unit of life energy.
"""
function walk!( natalia::NarrativeAgent, ne)
	wiggle!( natalia)
	move_agent!( natalia, ne, rand()*ne.search_speed)
	natalia.energy -= ne.living_cost
end

#--------------------------------------------------------------------------------------------------
"""
	explore!( natalia, ne)

Agent explores developmental possibilities to gain energy and so live longer.
"""
function explore!( natalia::NarrativeAgent, ne)
	if ne.narrative && !natalia.fixed
		# Narrative exploration:
		natalia.phenotype[natalia.p_genome] = rand(ne.alleles,count(natalia.p_genome))
		set_dissonance!( natalia, ne)
	end

	# Compare natalia's dissonance with that of an impressive nearby agent:
	nbrs = nearby_agents(natalia, ne, ne.interaction_distance)
	if !isempty(nbrs)
		nbrs = collect(nbrs)
		nbr = nbrs[findmax(ag->objective(ag,ne),nbrs)[2]]
		if objective(natalia,ne) > objective(nbr,ne)
			# Reward natalia for having locally lower dissonance:
			natalia.energy += ne.benefit * ne.living_cost
			natalia.fixed = true
		end
	end
end

#--------------------------------------------------------------------------------------------------
"""
	reproduce!(natalia, ne)

Mummy and Daddy agent give birth to two babies using mutation und crossover.
"""
function reproduce!(mummy::NarrativeAgent, ne)
	if nagents(ne) < ne.max_pop && mummy.energy > ne.birth_cost
		# I've got enough energy and there's sufficient space to have kids:
		mating_pool = nearby_agents(mummy, ne, ne.interaction_distance)
		if !isempty(mating_pool)
			mating_pool = collect(mating_pool)
			daddy = mating_pool[findmax(ag->objective(ag,ne),mating_pool)[2]]
			if daddy.energy > ne.birth_cost
				# Create progeny ...
				for (v_genome, p_genome) in crossover( mummy, daddy)
					θ = 2π * rand()
					baby = add_agent!( ne, (cos(θ),sin(θ)), v_genome, p_genome, copy(v_genome),
						ne.life_energy, 1.0, 1.0, 1.0, :black, false
					)
					mutate!( baby, ne)
					set_dissonance!( baby, ne)
				end

				# ... then deplete parents' energy due to birth:
				mummy.energy -= ne.birth_cost
				daddy.energy -= ne.birth_cost
			end
		end
	end
end

#--------------------------------------------------------------------------------------------------
"""
	crossover( mummy, daddy)

Perform crossover of value- and plasticity-genomes of two agents mummy and daddy.
"""
function crossover( mummy::NarrativeAgent, daddy::NarrativeAgent)
	# Select crossover points for both genomes:
	xpt_v, xpt_p = rand(1:length(mummy.v_genome)-1, 2)

	# Generate and return progeny's two genomes:
	[(
		vcat( daddy.v_genome[1:xpt_v], mummy.v_genome[xpt_v+1:end]),
		vcat( daddy.p_genome[1:xpt_p], mummy.p_genome[xpt_p+1:end])
	), (
		vcat( mummy.v_genome[1:xpt_v], daddy.v_genome[xpt_v+1:end]),
		vcat( mummy.p_genome[1:xpt_p], daddy.p_genome[xpt_p+1:end])
	)]
end

#--------------------------------------------------------------------------------------------------
"""
	mutate!( natalia, ne)

Perform mutation of value- and plasticity-genomes of the NarrativeAgent natalia.
"""
function mutate!( natalia::NarrativeAgent, ne)
	len = length(natalia.v_genome)

	# Mutate natalia's two genomes in place:
	mutated_loci = rand(len) .< ne.mu_rate
	natalia.v_genome[mutated_loci] = rand( ne.alleles, count(mutated_loci))
	mutated_loci = rand(len) .< ne.mu_rate
	natalia.p_genome[mutated_loci] = rand( Bool, count(mutated_loci))

	return natalia
end

#-----------------------------------------------------------------------------------------
"""
	set_dissonance!( natalia, ne)

Set natalia's dissonance and colour to reflect her dissonance with respect to the ne model's
target. Dissonance is measured in terms of the Hamming distance of natalia's genome from the
target, per available allele, per genome locus, and is guaranteed to lie in the range [0,1].
"""
function set_dissonance!( natalia::NarrativeAgent, ne)
	# Calculate dissonance based on phenotype and genotype:
	diss = sum(abs.(natalia.phenotype .- ne.target))/ne.problem_size
	natalia.plasticity = count(natalia.p_genome) / length(natalia.p_genome)
	if !ne.unimodal_fitness
		# Use bimodal dissonance function:
		diss = diss + 0.2sin(4pi*diss)
	end
	natalia.dissonance = min(diss, 1.0)
	if natalia.dissonance < natalia.min_diss
		natalia.min_diss = natalia.dissonance
	end

	# Colour natalia to indicate dissonance level:
	if natalia.dissonance > 0.5
		natalia.colour = :black
	elseif natalia.dissonance > 0.2
		natalia.colour = :blue
	elseif natalia.dissonance > 0.02
		natalia.colour = :red
	elseif natalia.dissonance > 0.002
		natalia.colour = :orange
	else
		natalia.colour = :yellow
	end
end

#--------------------------------------------------------------------------------------------------
"""
	objective( natalia, ne)

Calculate objective function of natalia within the given NarraticeExploration model.
"""
function objective( natalia::NarrativeAgent, ne)
	spiky(1-natalia.dissonance, ne.difficulty)
end

#--------------------------------------------------------------------------------------------------
"""
	create_problem(nalpha::Int, ntarget::Int)

Create a target string of length ntarget, constructed from symbols in an alphabet of nalpha
available alleles for the NarrativeExploration simulation. If nalpha < 2, the target is
constructed from the first ntarget characters of the AgentTools quotation, and the alphabet is the
set of unique characters in that target. If ntarget is greater than the length of the quotation,
the entire quotation is used as the target and the alphabet is a complete set of text characters.
If nalpha >= 2, the first nalpha letters of the English alphabet are used as the alphabet, and the
benchtest target is a string of length ntarget consisting of the first allele in that alphabet. For
nalpha == 2 and ntarget == 20, this reproduces the original Hinton and Nowlan (1987) problem.
"""
function create_problem(nalpha::Int=2, ntarget::Int=20)
	if nalpha < 2
		if ntarget < length(quotation)
			target = quotation[1:ntarget]
			return (collect(Set(target)), target)
		else
			return (quot_alphabet, quotation)
		end
	else
		alphabet = benchtest_alphabet[1:min(nalpha, length(benchtest_alphabet))]
		return (alphabet, fill(alphabet[1], ntarget))
	end
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Demonstriere den NarrativeExploration Versuch.
"""
function demo()
	# Parameters for interactive exploration:
	params = Dict(
		:narrative			=> [false, true],
		:mu_rate			=> 0:0.0005:0.005,
		:benefit			=> 1.0:0.1:3.0,
		:unimodal_fitness	=> (false, true),
		:difficulty			=> 0:0.1:1,
		:allele_count		=> 1:10,
		:target_length		=> 5:length(quotation),
		:life_energy		=> 1:10,
	)
	plotkwargs = (
		agent_color=(natalia->natalia.colour),
		agent_size=20,
		dt=1:200,
		mdata = [(m->m.mean_diss), (m->m.mean_plasticity), (m->m[m.min_dissID].min_diss)],
		mlabels = ["Mean dissonance", "Mean plasticity", "Minimum dissonance"],
	)

	playground,abmplt = abmplayground( narrative_exploration; params, plotkwargs...)

	# Dynamic observation of current lowest dissonance pheotype:
	best_v_genome = lift( wld -> "v: "*String(wld[wld.min_dissID].v_genome), abmplt.model)
	best_p_genome = lift((wld -> "p: "*join(Int.(wld[wld.min_dissID].p_genome))), abmplt.model)

	# Display as text in the plot:
	text!( 10,1.4, text=best_p_genome, color=:blue,
		font = "Courier New Bold", fontsize = 16, align = (:left,:top))
	text!( 10,1.3, text=best_v_genome, color=:red,
		font = "Courier New Bold", fontsize = 16, align = (:left,:top))

	playground
end

end     # ... of module NarrativeExploration 