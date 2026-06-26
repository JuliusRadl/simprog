#========================================================================================#
"""
	TextAnalysis

Module TextAnalysis: This is a collection of methods for analysing text.

Author: Niall Palfreyman, 01/08/2023
"""
module TextAnalysis

# Externally callable symbols of the TextAnalysis module:
export splitwords, entry_counts, ngrams, ngram_counts, completion_cache

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	splitwords( text)

Generate a list of words in the text, assuming the usual word-boundary separators.
"""
function splitwords( text)
	cleantext = replace( text, r"\s+" => " ")	# Clean up whitespace
	words = split( cleantext, r"(\s|\b)")		# Split words on whitespace or other word boundary

	words
end

"""
	entry_counts( list::Vector)

Generate a dictionary containing all entries in the entrylist, accompanied by their frequency
in the entrylist.
Arguments: A Vector
Values: A dictionary where the keys are unique vector elements
and the values are how often the element appears 
in the original vector
"""
function entry_counts(list::Vector)
        counts = Dict()
        for entry in list
            if haskey(counts,entry)
                counts[entry] += 1
            else
                counts[entry] = 1
            end
        end
        counts
    end

"""
	ngrams( wordlist::Vector, n::Int)

Construct a list of all word-sequences of length n contained in the given wordlist.
"""
function ngrams( wordlist::Vector, n::Int)
	startingPositions = 1:length(wordlist)-n+1
	map(startingPositions) do startingPosition
		wordlist[startingPosition:startingPosition+n-1]
	end
end

"""
	ngram_counts( wordlist::Vector, n::Int)

Construct a Dictionary of counts of all word-sequences of length n in the given wordlist.
"""
function ngram_counts(wordlist::Vector, n::Int)
	entry_counts(ngrams(wordlist, n))
end

"""
	completion_cache( grams)

Construct a Dictionary of all word-sequences of length n-1, and assign to each such key a ROW 
vector containing words that complete this to an n-gram in the given collection of grams.
"""
function completion_cache(grams)
        cache = Dict()
        for ngram in grams
            if haskey(cache,ngram[1:end-1])
				# Append last word to cache row vector
                cache[ngram[1:end-1]] = [cache[ngram[1:end-1]] ngram[end]]
            else
                cache[ngram[1:end-1]] = [ngram[end]]
            end
        end

        cache
    end

#-----------------------------------------------------------------------------------------
"""
	demo()

Your module should have a demo() method at the end to show users how to use your module.
Notice that I haven't exported the demo() method - users should call it explicitly like this:
	TextAnalysis.demo()
"""
function demo()
	println("\n============ Demonstrate TextAnalysis: ===============")
	println("First create a sampleText:")
	sampleText = "A lazy   brown fox  trips,           over the    lazy brown dog."
	display( sampleText)
	println()

	word_list = splitwords(sampleText)
	println("Now divide up the sampleText into individual words:")
	display( word_list)
	println()

	println("Here is a dictionary of the word frequencies:")
	display( entry_counts(word_list))
	println()

	println("Now lets create the corresponding trigrams")
	sampleTextTrigrams = ngrams(word_list, 3)
	display(sampleTextTrigrams)
	println()

	println("How often do each of these trigrams appear?")
	sampleTextTrigramCounts = ngram_counts(word_list, 3)
	display(sampleTextTrigramCounts)
	println()

	println("Let's create a dictionary in which the keys are the words of the ngram except",
	"for the last one, which are the values:")
	cache = completion_cache(sampleTextTrigrams)
	display(cache)
	println()

end

end