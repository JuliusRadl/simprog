#========================================================================================#
"""
	Sumprods

Module Sumprods: Rewrite this module IN YOUR OWN CODE in good julian style to demonstrate and
explain in detail which of the following methods is the BEST implementation of sumprod().

Source: https://discourse.julialang.org/t/julian-way-to-write-this-code/119348

Author: Niall Palfreyman, 16/09/2024
"""
module Sumprods

using BenchmarkTools

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"sumprod1(a,b,c): Sum the products of corresponding elements of three vectors a, b, c."
function sumprod1(a,b,c)
	sum(i -> a[i] * b[i] * c[i], 1:length(a))
end

"sumprod2(a,b,c): Sum the products of corresponding elements of three vectors a, b, c."
function sumprod2(a,b,c)
	sum(a .* b .* c)
end

"sumprod3(a,b,c): Sum the products of corresponding elements of three vectors a, b, c."
function sumprod3(a,b,c)
	a .*= b
	a .*= c
	sum(a)
end

"sumprod4(a,b,c): Sum the products of corresponding elements of three vectors a, b, c."
sumprod4(a,b,c) = sum(splat(*), zip(a,b,c))

"""
Demonstrate usage of the four sumprod() implementations.
"""
function demo()
	p,q,r = [rand(4) for _ in 1:3]

	display( @benchmark sumprod1($p,$q,$r))
	display( @benchmark sumprod2($p,$q,$r))
	display( @benchmark sumprod3($p,$q,$r))
	display( @benchmark sumprod4($p,$q,$r))
end

end