#========================================================================================#
"""
	Permutator

Utility Functions for generating and applying permutation vectors.

Author: Julius Radl
"""

module Permutator

export get_permutations, apply_permutations

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
    rec_permutations(elements, permutation, permutations, elements_left)

Helper function for find_permutations().
Recursively find all permutations of elements in an array.
"""
function rec_permutations(num_list, permutation, perm_vectors, nums_left)
    if nums_left <= 0
            push!(perm_vectors, copy(permutation))
        return
    else
        for i in 1:length(num_list)
            num = num_list[i]
            deleteat!(num_list, i)
            push!(permutation, num)

            rec_permutations(num_list, permutation, perm_vectors, nums_left-1)

            pop!(permutation)
            insert!(num_list, i, num)
        end
    end
end

#-----------------------------------------------------------------------------------------
"""
    find_permutations(elements::Vector{Any})

Generate all permutation vectors for num_elements elements.
    
We don't need a method of this function for other data
types, as a list of elements of any type can be permuted
using a permutation vector containg integers 
(see apply_permutations()).
"""
function get_permutations(num_elements::Int64)
    perm_vectors = Vector{Vector{Int64}}()
    permutation = Vector{Int64}()

    num_list = collect(1:num_elements)

    rec_permutations(num_list, permutation, perm_vectors, num_elements)
    return perm_vectors
end

#-----------------------------------------------------------------------------------------
"""
    apply_permutations(elements, perm_vectors)

Get all permutations of the vector elements using a
vector perm_vectors containing permutation vectors.
"""
function apply_permutations(elements, perm_vectors::Vector{Vector{Int64}})
    # Type the output vector according to the input
    permutations = Vector{typeof(elements)}()

    for i in 1:length(perm_vectors)
        perm = elements[perm_vectors[i]]
        push!(permutations, perm)
    end
    
    return permutations
end

#-----------------------------------------------------------------------------------------
function demo()
    println("Let's demonstrate how to use Permutators.jl!\n")

    println("First we need a list of elements. It may contain elements ",
    "of any type. I'll choose Strings:")
    elements = ["apples", "oranges", "bananas"]
    display(elements)
    println()

    println("Now let's get all permutation vectors for a list of this length:")
    len = length(elements)
    println("Length of our list: ", len)
    perm_vecs = get_permutations(len)
    display(perm_vecs)
    println()

    println("We are ready to apply the permutation vectors to our list:")
    permutations = apply_permutations(elements, perm_vecs)
    display(permutations)
    println()
end

end # End of Permutations.jl

