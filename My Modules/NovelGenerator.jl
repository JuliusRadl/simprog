#========================================================================================#
"""
	NovelGenerator

Module NovelGenerator: This is a collection of methods for analysing and generating text.

Author: Julius Radl, 08/11/2025
Based on Code by Niall Palfreyman
"""
module NovelGenerator

# Imports
import Dates
import Printf

# Externally callable symbols of the NovelGenerator module:
export splitwords, entry_counts, ngrams, ngram_counts, completion_cache, getText, trimText

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	getText(txtPath::String)

Get all text from a .txt-File as a String.
"""
function getText(txtPath::String)
	textData = open(txtPath, "r")
	rawText = read(textData, String)
	close(textData)

	rawText
end

"""
	trimText(beginningPhrase::String, endingPhrase:: String)

Trim a text String to the desired content
"""
function trimText(text::String, beginningPhrase::String, endingPhrase:: String)
	beginning_index = findfirst(beginningPhrase, text)[1]
	ending_index = findlast(endingPhrase, text)[end]
	trimmedText = text[beginning_index : ending_index]

	trimmedText
end

"""
	splitwords(text::String)

Generate a list of words in the text, assuming the usual word-boundary separators.
"""
function splitwords(text::String)
	cleantext = replace(text, r"\s+" => " ")	# Clean up whitespace
	words = split(cleantext, r"(\s|\b)")		# Split words on whitespace or other word boundary

	words
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
	completion_cache(ngrams)

Construct a Dictionary of all word-sequences of length n-1, and assign to each such key a ROW 
vector containing words that complete this to an n-gram in the given collection of grams.
"""
function completion_cache(ngrams)
        cache = Dict()
        for ngram in ngrams
            if haskey(cache,ngram[1:end-1])
				# Append last word to cache row vector
                cache[ngram[1:end-1]] = [cache[ngram[1:end-1]] ngram[end]]
            else
                cache[ngram[1:end-1]] = [ngram[end]]
            end
        end

        cache
    end

"""
	write_novel(source_text::String, num_words::Int, n=3)

Generate new text in the style of the input text using a language model.
Save it to a file.
"""
function write_novel(sourceText::String, num_words::Int, n=3)
	# Set up the language model
	wordList = splitwords(sourceText)
	circularWordList = [wordList..., wordList[1:n-1]...]
	nGrams = ngrams(circularWordList, n)
	cache = completion_cache(nGrams)

	# Start the novel with a capital letter
	novelWordList = rand(nGrams)
	while !occursin(r"[\p{Lu}“]", novelWordList[1])
		novelWordList = rand(nGrams)
	end

	# Fill the novel and end it with a sentece-ending special character
	while !occursin(r"[!?.”]", novelWordList[end]) || length(novelWordList) < num_words
		previous_words = novelWordList[end-(n-2):end]
		possible_words = cache[previous_words]
        new_word = rand(possible_words)
        push!(novelWordList, new_word)
	end

	novelText = join_wordlist(novelWordList)

	# Construct a header
	copyright = 
	title = join(split(novelText)[1:3], " ")
	header = (
		"#----------------------------------------------#\n" *
		Printf.@sprintf("|%30s                |\n", title) *
		"#----------------------------------------------#\n\n"
	)

	# Finish the novel
	novelText = (getCopyright() * "\n\n" * header * novelText)
	
	novelText
end

"""
	join_wordlist(wordlist::Vector{String})

Joins the words in a wordlist smoothly into a semi-coherent text.
"""
function join_wordlist(wordlist::Vector)

	joined_wordlist = join(wordlist, " ")

	joined_wordlist = replace(joined_wordlist, " ," => ",")
	joined_wordlist = replace(joined_wordlist, " ;" => ";")
	joined_wordlist = replace(joined_wordlist, " ”" => "”")
	joined_wordlist = replace(joined_wordlist, " ." => ".")
	joined_wordlist = replace(joined_wordlist, " !" => "!")
	joined_wordlist = replace(joined_wordlist, " ?" => "?")
	joined_wordlist = replace(joined_wordlist, "“ " => "“")
	joined_wordlist = replace(joined_wordlist, " ’ " => "’")
	

	joined_wordlist
end

"""
	getCopyright(author::String)

Creates a copyright string.
"""
function getCopyright()
	
	date = Dates.today()
	dateString = Dates.format(date, "yyyy-mm-dd")
	copyright = ("--- © " * dateString * " ---")
	
	copyright
end
	
#-----------------------------------------------------------------------------------------
"""
	demo()

Your module should have a demo() method at the end to show users how to use your module.
Notice that I haven't exported the demo() method - users should call it explicitly like this:
	NovelGenerator.demo()
"""
function demo()
	println("\n============ Demonstrate NovelGenerator: ===============\n")
	println("Let's generate a new text based on a sample text! For this to work " *
	"you should have a file called pandp_text.txt in the same directory as NovelGenerator.jl!")
	
	println("\nTime to create the path String:")
	dirPath = @__DIR__
	filePath = joinpath(dirPath, "pandp_text.txt")
	display(filePath)
	
	println("\nNow we can extract the text:")
	rawText = getText(filePath)
	display(rawText)

	println("\nLet's trim the text: Pride and Prejudice's actual story " *
	"starts with the words 'It is a truth' and ends with 'uniting them.':")
	trimmedText = trimText(rawText, "It is a truth", "uniting them.")
	display(trimmedText)

	println("\nNow we are ready to generate a new text based on Pride and Prejudice. " *
	"Let's make it 200 Words long:\n")
	newText = write_novel(trimmedText, 200, 3)
	print(newText)

end

end