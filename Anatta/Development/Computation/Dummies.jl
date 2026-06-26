#========================================================================================#
"""
	Dummies

Module Dummies: This is a dummy module. It contains a data type Dummy, plus a list of methods
that users can use to work with your data type. Use this file as a template for creating your
own julia programs.

Author: Niall Palfreyman, 05/07/2023
"""
module Dummies

# Externally callable symbols of the Dummies module:
export Dummy, greeting

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Dummy

This is a dummy data type.
"""
struct Dummy
	name::String			# The name of the Dummy
	age::Int				# The age of the Dummy

	function Dummy( name::String, age)
		int_age = round(Int,age)
		return new(name,int_age)
	end
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	greeting( dummy)

Create a greeting from the Dummy. This is just an example method to show you how to write your own
methods for manipulating the Dummy data type defined above.
"""
function greeting( dummy::Dummy)
	string( "Hello, I am a Dummy named ", dummy.name, ". I am ", dummy.age, " years old! :)")
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Your module should have a demo() method at the end to show users how to use your module.
Notice that I haven't exported the demo() method - users should call it explicitly like this:
	Dummies.demo()
"""
function demo()
	println("\n============ Demonstrate Dummies: ===============")
	println("First create a Dummy:")
	dimmy = Dummy( "Dimmy", 3.1415)
	display( dimmy)
	println()

	println("Now display the Dummy's greeting:")
	println( greeting(dimmy))
	println()
end

end