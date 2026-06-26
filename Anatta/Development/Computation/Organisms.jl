#========================================================================================#
"""
	Organisms

Module Organisms: This is a Weasel module. It contains a data type Weasel, plus a list of methods
that users can use to work with your data type. Use this file as a template for creating your
own julia programs.

Author: Niall Palfreyman, 05/07/2023
"""
module Organisms

# Externally callable symbols of the Organisms module:
export Weasel, greeting, Organism, Animal, Rabbit, meet, encounter

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Organism

A general living being
"""
abstract type Organism end

#-----------------------------------------------------------------------------------------
"""
	Animal

An animate Organism
"""
abstract type Animal <: Organism end

#-----------------------------------------------------------------------------------------
"""
	Weasel

This is a Weasel data type.
"""
struct Weasel <: Animal
	name::String			# The name of the Weasel
	age::Int				# The age of the Weasel

	function Weasel( name::String, age)
		int_age = round(Int,age)
		return new(name,int_age)
	end
end

#-----------------------------------------------------------------------------------------
"""

	Rabbit

This is a Rabbit data type.
"""

mutable struct Rabbit <: Animal
	name::String
    age::Integer
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	greeting(Animal)

Create a greeting from the Organism. This is just an example method to show you how to write your own
methods for manipulating the Weasel data type defined above.
"""
function greeting(animal::Animal)
	string("Hello, I am a ", string(nameof(typeof(animal))), " named ", animal.name, ". I am ", animal.age, " years old! :)")
end

#-----------------------------------------------------------------------------------------
"""
	meet()

Return a string based on the (sub-)types of two organisms,
indicating how these organisms would interact during a meeting.
"""
meet( meeter::Weasel, meetee::Rabbit) = "attacks"
meet( meeter::Weasel, meetee::Weasel) = "challenges"
meet( meeter::Rabbit, meetee::Rabbit) = "sniffs"
meet( meeter::Rabbit, meetee::Weasel) = "hides"
meet( meeter::Organism, meetee::Organism) = "ignores"

#----------------------------------------------------------------------------------------
"""
	encounter(meeter<:Organism, meetee<:Organism)

Print a line describing a meeting between two organisms using the meet-Function,
which has multiple methods depending on the organisms' subtypes
"""

function encounter(meeter::Organism, meetee::Organism)
	println(meeter.name, " meets ", meetee.name, " and ", meet(meeter, meetee), ".")
end
	
#-----------------------------------------------------------------------------------------
"""
	demo()

Your module should have a demo() method at the end to show users how to use your module.
Notice that I haven't exported the demo() method - users should call it explicitly like this:
	Organisms.demo()
"""
function demo()
	println("\n============ Demonstrate Organisms: ===============")
	println("First create a Weasel:")
	wendy = Weasel( "Wendy", 3.1415)
	display( wendy)
	println()

	println("Now display the Weasel's greeting:")
	println( greeting(wendy))
	println()

	println("Now let's create two Rabbits")
	rabia = Rabbit("Rabia", 5)
	display(rabia)
	println("Now display the first Rabbit's greeting:")
	println(greeting(rabia))
	println("Let's create the second Rabbit:")
	robby = Rabbit("Robby", 34)
	display(robby)
	println("Now display the second Rabbit's greeting:")
	println(greeting(robby))
	println()

	println("Now let's have the rabbits meet:")
	encounter(rabia, robby)
	println("But what if a rabbit meets the weasel?")
	encounter(rabia, wendy)
	println()

	println("Now isnt' that just neat?!")
	println()

	println("We can, of course, create a new type of organism, a tree:")
	codeString = "struct Tree <: Organism; name::String; age::Int; end"
	println(codeString)
	eval(Meta.parse(codeString))
	println()

	println("Why don't we grow a tree?")
	tongtong = Tree("Tongtong", 300)
	display(tongtong)
	println("Naturally, trees can't greet.")
	println()

	println("Let's have a rabbit meet a tree:")
	encounter(robby, tongtong)
end

end