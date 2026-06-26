#========================================================================================#
"""
	Utilities

A library of useful functions.

Author: Niall Palfreyman, 06/09/2022.
"""
module Utilities

export mean, nestlist, normed, sample, std

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	eratosthenes_bad( N)

Generate a list of prime numbers up to a user specified maximum N. This implementation ports
directly from Java to Julia, and is therefore unreasonably low-level.
"""
function eratosthenes_bad( N)
	if N >= 2 # the only valid case
		# declarations:
		f = Bool[]
		i = 0
		# initialise array to true
		for i in 1:N
			push!(f, true)
		end

		# get rid of known non-primes
		f[1] = false

		# sieve
		j = 0
		for i in 2:round(Int, sqrt(N)+1)
			if f[i] # if i is uncrossed, cross its multiples
				for j in 2*i:i:N
					f[j] = false # multiple is not a prime
				end
			end
		end

		# how many primes are there
		count = 0
		for i in 1:N
			if f[i]
				count += 1
			end
		end

		primes = zeros(Int, count)

		# move the primes into the result
		j = 1
		for i in 1:N
			if f[i] # if prime
				primes[j] = i
				j += 1
			end
		end
		return primes # return the primes
	else # if N < 2
		return Int[] # return null array if bad imput
	end
end

#-----------------------------------------------------------------------------------------
"""
	eratosthenes_good( N)

Generate a list of prime numbers up to a user specified maximum N. This implementation ports
directly from Java to Julia, and is therefore unreasonably low-level.

This implementation is identical with that of eratosthenes_bad() - it is your job to rewrite it
into GOOD julian code.
"""
function eratosthenes_good( N)
	if N >= 2 # the only valid case
		# declarations:
		f = Bool[]
		i = 0
		# initialise array to true
		for i in 1:N
			push!(f, true)
		end

		# get rid of known non-primes
		f[1] = false

		# sieve
		j = 0
		for i in 2:round(Int, sqrt(N)+1)
			if f[i] # if i is uncrossed, cross its multiples
				for j in 2*i:i:N
					f[j] = false # multiple is not a prime
				end
			end
		end

		# how many primes are there
		count = 0
		for i in 1:N
			if f[i]
				count += 1
			end
		end

		primes = zeros(Int, count)

		# move the primes into the result
		j = 1
		for i in 1:N
			if f[i] # if prime
				primes[j] = i
				j += 1
			end
		end
		return primes # return the primes
	else # if N < 2
		return Int[] # return null array if bad imput
	end
end

#-----------------------------------------------------------------------------------------
"""
	senehtsotare(N)

Generate a list of prime numbers up to a user specified maximum N. This is George Dateris'
'good' implementation, which incorporates the following changes:
	* Replace initial `if` clause by an early return statement;
	* Crossing out integers is implemented within its own function;
	* Uses Julia's efficient search function for true elements in a boolean array;
	* Uses Julia's 1-based array indexing of `isprime` for efficient crossing-out.
"""
function senehtsotare(N::Int)
	N < 2 && return Int[]				# By definition, there are no primes less than 2

	isprime = trues(N)					# Number `N` is prime if `isprime[N] == true`
	isprime[1] = false					# By definition, 1 is not prime
	crossout_prime_multiples!(isprime)

	findall(isprime)
end

#-----------------------------------------------------------------------------------------
"""
	crossout_prime_multiples!( isprime::AbstractVector{Bool})

For all primes in `isprime` (elements that are `true`), set all their multiples to `false`.
Assumes `isprime` starts counting from 1.
"""
function crossout_prime_multiples!( isprime::AbstractVector{Bool})
	N = length(isprime)
	for i in 2:round(Int, sqrt(N)+1)
		if isprime[i]
			for j in 2i:i:N
				isprime[j] = false
			end
		end
	end
end

#-----------------------------------------------------------------------------------------
"""
	nestlist( f, x0, n)

Apply the function f repeatedly to the initial value x0, creating a list of n values generated
in this way. This is useful for simulating the development of initial data over time when repeatedly
acted upon by a discrete rule.
"""
function nestlist( f::Function, x0, n::Integer)
	if n ≤ 0
		# Base case. Return the initial value x0:
		[x0]
	else
		# Recurrence case. Fill a list with n applications of f to the initial value x0:
		list = Vector{typeof(x0)}(undef, n+1)			# List has correct size and type
		list[1] = x0									# Initial value
		for i in 1:n
			# Enter i-th iterated value:
			list[i+1] = f(list[i])
		end
		list											# Return filled list
	end
end

#-----------------------------------------------------------------------------------------
"""
	normed( v::Vector{Float64})

Return a vector of components in the same proportions as the elements of the input vector v, but
normalised to sum to 1.0. Useful for converting a vector of frequencies into a probability.
"""
normed( v::Vector{Float64}) = v/sum(v)

#-----------------------------------------------------------------------------------------
"""
	mean( v::Vector{Float64})

Return the mean of the elements of the input vector v.
"""
mean( v::Vector{Float64}) = sum(v)/length(v)

#-----------------------------------------------------------------------------------------
"""
	std( v::Vector{Float64})

Return the standard deviation of the elements of the input vector v.
"""
std( v::Vector{Float64}) = sqrt(sum((v.-mean(v)).^2)/max(1,length(v)-1))

#-----------------------------------------------------------------------------------------
"""
	sample( x::Vector{Float64}, n=length(x), T=0.0)

Create a vector of n sampled types from a population whose types have population frequencies
given by the vector x. Sampling is implemented using roulette-wheel selection based on a
roulette wheel whose pockets are biased according to the frequencies.

T is a "Temperature" parameter that determines the extent to which sampling depends upon the
frequency distribution. When T is high (T>=1.0), the sampling is (rigidly) uniformly distributed;
when T is low (T<=0.0), sampling becomes increasingly determined by the frequency distribution
vector x. This implementation also uses sigma-scaling to eliminate from sampling all types whose
frequency is smaller than 1/(2T+0.1) standard deviations below the frequency mean.
"""
function sample( x::Vector{Float64}, n=length(x); T=0.0)
	# Preprocess frequencies:
	T = min(1.0,max(0.0,T))								# Enforce 0.0 <= T <= 1.0
	x = T*mean(x) .+ (1-T)*x							# Smooth frequencies towards average
	x = x / sum(x)										# Normaise the resulting frequencies
	avg = mean(x)										# Average population frequency
	stdev = std(x)										# Spread of frequencies from average

	# Create fitness values:
	if stdev == 0										# Singular case: All frequencies were
		fitnesses = ones(Float64,length(x))				# equal to the average!
	else
		fitnesses = 1.0 .+ (x.-avg)/((0.1+2T)*stdev)	# Frequencies smaller than 1/(2T+0.1) std
		fitnesses[fitnesses .<= 0] .= 0.0				# deviations below average are unviable
	end
	
	# Create fitness-biased roulette wheel:
	wheel = cumsum( fitnesses)/sum(fitnesses)			# Roulette-wheel of fitness-biased pockets
	throws = (rand() .+ (1:n)/n) .% 1.0					# Stochastic Universal Sampling: n ball-throws
														# distributed uniformly over interval (0,1)
	# Create vector of n sampled population types:
	samples = zeros(Integer,n)
	for (throw,sample) in enumerate(throws)				# Ball-throw delivers one random sample fitness
		for (pocket,fitness) in enumerate(wheel)		# Wheel pocket contains fitness of one type
			if sample <= fitness						# When the sample surpasses type's fitness ...
				samples[throw] = pocket					# ... use this pocket as the sampled type
				break									# ... then skip on to the next throw.
			end
		end
	end

	(n==1) ? samples[1] : samples						# Strip brackets if we only want one sample
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Use-case demonstrations of the method sample().
"""
function demo( temperature=0.0)
	x = 1:10.0
	x = collect(x/sum(x))

	println( "Sampling the frequencies")
	println( "    $(round.(x,digits=2))")
	println( "generates the following sampled types:")
	println( "    ", sample(x,T=temperature))
	println()
	println( "Sampling only 7 types generates:")
	println( "    ", sample(x,7,T=temperature))
end

end		# of Utilities