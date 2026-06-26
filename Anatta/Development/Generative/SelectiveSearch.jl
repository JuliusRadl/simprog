#========================================================================================#
"""
	Selective Search

A simulation to investigate the selective searching of genetic structures. SelectiveUnits do not
move, and do not really interact with each other. Rather, they are passive storage units in a
single Genetic Algorithm machine that manipulates their genomes based on how well they minimise
the dissonance between their genome and a target string. The SelectiveUnits' genomes are character
strings, and their dissonance is the mean Hamming distance per character of their genome from the
target quotation string from AgentTools. The simulation aims to evolve the population of
SelectiveUnits over successive generations towards the target string.

Author: Niall Palfreyman, March 2025.
"""
module SelectiveSearch

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, Random, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	SelectiveUnit

SelectiveUnits sit in one place, unmoving, and are selected for reproduction depending on how well
they minimise a particular objective function: dissonance.
"""
@agent struct SelectiveUnit(ContinuousAgent{2,Float64})
	genome::Vector{Char}			# SelectiveUnit's character genome
	dissonance::Float64				# Dissonance of the SelectiveUnit with its environment
	colour::Symbol					# Colour indicating SelectiveUnit's dissonance status
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	selective_search()

Initialise the Selective Search simulation.
"""
function selective_search(;
	crossover = true,
	mu_rate = 1e-3,
	elitism = 0.0,
	unimodal_fitness = true,
)
	extent = (40,40)									# Extent of the world
	nPop = round(Int,0.3*prod(extent))					# Population size
	properties = Dict(
		:crossover			=> crossover,				# Are we using crossover?
		:mu_rate			=> mu_rate,					# Probability of mutating a locus
		:elitism			=> elitism,					# How much dissonance are we eliminating?
		:unimodal_fitness	=> unimodal_fitness,		# Use unimodal fitness landscape?
		:nPop				=> nPop,					# (Constant) population size
		:roulette			=> cumsum(ones(nPop)/nPop),	# Roulette wheel selector
		:target				=> quotation,				# Target string to be evolved
		:alphabet			=> quot_alphabet,			# Collection of available alleles
		:mean				=> 1.0,						# Mean dissonance of population
		:minID				=> 1,						# ID of SelectiveUnit with minimum dissonance
	)

	ss = StandardABM( SelectiveUnit, ContinuousSpace(extent);
		model_step!, properties
	)
	genome_length = length(ss.target)
	for _ in 1:nPop
		susie = add_agent!( ss, (1,1),					# Ur-SelectiveUnit has ...
			rand( ss.alphabet, genome_length),			# Genome (candidate Vector{Char})
			1.0, :white									# Provisional dissonance and fitness type
		)
		set_dissonance!( susie, ss)						# Set susie's dissonance/colour
	end

	ss
end

#-----------------------------------------------------------------------------------------
"""
	model_step!( ss)

Perform dissonance-based selection on the SelectiveUnit population to decide who gets children.
"""
function model_step!( ss::ABM)
	allSUs = collect(allagents(ss))						# List of all SelectiveUnits
	dissonances = (ag->ag.dissonance).(allSUs)			# List of their dissonances
	ss.minID = allSUs[findmin(dissonances)[2]].id		# Record SelectiveUnit with minimum dissonance
	parents = select_parents( dissonances, ss)			# Select parents based on dissonance

	#-------------------------------------------------------------------------------------
	# Recombine all parent genomes to create next generation children genomes with mutation:
	nMatings = length(parents) ÷ 2						# N parents support N ÷ 2 matings
	genomes = (ag->ag.genome).(allSUs)					# List of all genomes
	gene_length = length(genomes[1])					# Length of each genome	
	mummies = genomes[parents[1:nMatings]]				# 1st half of parents are Mummies ...
	daddies = genomes[parents[nMatings+1:end]]			# 2nd half are Daddies
	children = similar(genomes)							# Preallocate Vector of child genomes
	
	for m in 1:nMatings
		xpt = rand(1:gene_length-1)						# Crossover point
		sally = deepcopy.(mummies[m])					# Sally is mummy-based child ...
		sally[xpt+1:end] = daddies[m][xpt+1:end]		# who also takes after daddy.
		billy = deepcopy.(daddies[m])					# Billy is daddy-based child ...
		billy[xpt+1:end] = mummies[m][xpt+1:end]		# who also takes after mummy.

		# Mutate both children of the mating:
		if ss.mu_rate > 0.0
			mutation_loci = rand(gene_length).<ss.mu_rate
			sally[mutation_loci] = rand(ss.alphabet,count(mutation_loci))
			mutation_loci = rand(gene_length).<ss.mu_rate
			billy[mutation_loci] = rand(ss.alphabet,count(mutation_loci))
		end

		# Store sally and billy in children genomes:
		children[m] = sally
		children[m+nMatings] = billy
	end

	#-------------------------------------------------------------------------------------
	# Update all SelectiveUnits with their new children's genomes:
	for (i,susie) in enumerate(allSUs)
		susie.genome = children[i]
		set_dissonance!( susie, ss)
	end
end

#-----------------------------------------------------------------------------------------
"""
	select_parents( dissonances, ss)

Perform dissonance-based roulette-wheel selection on the provided list of dissonances to decide
which of them are most likely to be selected as parents for the next generation.
"""
function select_parents( dissonances, ss)
	ss.mean = mean(dissonances)
	sigma = std(dissonances)							# Standard deviation

	#-------------------------------------------------------------------------------------
	# Create roulette-wheel whose slot-widths reflect dissonance-based fitness of each parent:
	if sigma <= 1e-6
		# Singular case: all evaluations were approximately equal to the mean:
		fitnesses = ones(ss.nPop)
	else
		# Otherwise use sigma-scaling and exclude bottom elitism tail of population distribution:
		fitnesses = (ss.mean .- dissonances)/sigma .- ss.elitism
		fitnesses[fitnesses .<= 0] .= 0
	end
	roulette = cumsum( fitnesses)/sum(fitnesses)		# Roulette wheel of dissonance-biased slots
	throws = rem.( rand() .+ (1:ss.nPop)./ss.nPop, 1)	# Statistical Uniformly Sampled ball throws

	#-------------------------------------------------------------------------------------
	# Roulette-wheel selection - Create randomised list of fitness-selected parent genome indices
	# by successively throwing the roulette ball onto the wheel and recording the slot it lands on:
	parents = zeros(Integer,ss.nPop)					# List of fitness-selected parents
	for parent in 1:ss.nPop
		for slot in 1:ss.nPop
			if throws[parent] <= roulette[slot]
				parents[parent] = slot
				break
			end
		end
	end

	shuffle(parents)									# Randomise parent order
end

#-----------------------------------------------------------------------------------------
"""
	set_dissonance!( susie, ss)

Set SelectiveUnit's initial colour to reflect its dissonance with respect to the ss model's
target string. A SelectiveUnit's dissonance is the mean Hamming distance per character of its genome from
the target.
"""
function set_dissonance!( susie::SelectiveUnit, ss)
	if ss.unimodal_fitness
		susie.dissonance = dissonance( susie.genome, ss)	# Mean Manhattan difference per locus
	else
		susie.dissonance = bimodal_dissonance( susie.genome, ss)
	end

	# Colour SelectiveUnit to indicate dissonance level
	if susie.dissonance > 0.5
		susie.colour = :black
	elseif susie.dissonance > 0.2
		susie.colour = :blue
	elseif susie.dissonance > 0.02
		susie.colour = :red
	elseif susie.dissonance > 0.002
		susie.colour = :orange
	else
		susie.colour = :yellow
	end
end

#-----------------------------------------------------------------------------------------
"""
	dissonance( candidate::Vector{Char}, ss)

The dissonance of a candidate string with respect to a target string (of the same length as the
candidate) is its Hamming distance from the target, divided by the length of the target - so, the
mean Hamming distance per string character.
"""
function dissonance( candidate::Vector{Char}, ss)
	sum(abs.(candidate - ss.target))/(length(ss.target)*length(ss.alphabet)-1)
end

#-----------------------------------------------------------------------------------------
"""
	bimodal_dissonance( candidate::Vector{Char}, ss)

This dissonance function still has the same global minimum at 0 as dissonance() when the
candidate and target strings are equal to each other, but introduces a confusing local minimum.
"""
function bimodal_dissonance( candidate::Vector{Char}, ss)
	diss = dissonance(candidate, ss)
	diss + 0.2sin(4pi*diss)
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Demonstrate a simple genetic algorithm.
"""
function demo()
	params = Dict(
		:crossover => false:true,
		:mu_rate => 0:1e-4:0.01,
		:elitism => -2:0.1:0.2,
		:unimodal_fitness => (false, true),
	)
	plotkwargs = (
		agent_color=(SelectiveUnit->SelectiveUnit.colour),
		agent_size=20,
		mdata = [(m->m.mean), (m->m[m.minID].dissonance)],
		mlabels = ["Mean dissonance", "Minimum dissonance"],
	)

	playground,abmplt = abmplayground( selective_search; params, plotkwargs...)

	best_string = lift( (wld->String(wld[wld.minID].genome)), abmplt.model)
	text!( 10,1.0, text=best_string, color=:red, fontsize=16, align=(:left,:top))

	playground
end

end