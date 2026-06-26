#========================================================================================#
"""
	AgentTools

This module provides a collection of convenient utility functions to supplement the Agents package.

Author:  Niall Palfreyman, January 2025.
"""
module AgentTools

using Agents, GLMakie, GeometryBasics, Observables

import InteractiveUtils:@which

export abmplayground, multicoloured, dejong2, diffuse4, diffuse4!, diffuse8, diffuse8!,
		gradient, mean, norm, size, spectrum, std, turn!, left!, right!, wiggle!, valleys, wedge

#-----------------------------------------------------------------------------------------
# Module data:
#-----------------------------------------------------------------------------------------
"""
	spectrum

Spectrum of colours for displaying the id of agents.
"""
const spectrum = [:darkblue,:blue,:green,:violet,:crimson,:red,:orange,:yellow,:white]

#-----------------------------------------------------------------------------------------
"""
	wejj

Basic wedge shape for use in wedge() function.
"""
const wejj = [[1,0],[-0.5,0.5],[-0.5,-0.5]]

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	size( cs::ContinuousSpace)

Extend Base.size to ContinuousSpaces to allow heatmaps to work in abmplot.
"""
Base.size(cs::Agents.ContinuousSpace) = cs.dims

#-----------------------------------------------------------------------------------------
"""
	norm( arr)

Calculate the norm of the array or tuple arr.
"""
function norm( arr)
	sum(abs.(arr).^2)^0.5
end

#-----------------------------------------------------------------------------------------
"""
	mean(arr)

Calculate the mean value contained in the array or tuple arr.
"""
function mean( arr)
	sum(arr)/length(arr)
end

#-----------------------------------------------------------------------------------------
"""
	std(arr)

Calculate the standard deviation of the values contained in the array or tuple arr.
"""
function std( arr)
	avg = mean(arr)
	sqrt(mean((arr.-avg).^2))
end

#-----------------------------------------------------------------------------------------
"""
	multicoloured(agent)

Return different colours depending on the id of the agent.
"""
function multicoloured( agent::AbstractAgent)
	spectrum[1+agent.id%length(spectrum)]
end

#-----------------------------------------------------------------------------------------
"""
	wedge(agent::AbstractAgent, siz=1.0)

Draw a wedge-shaped marker pointing in the agent's facing (vel) direction.
"""
function wedge( agent::AbstractAgent, siz=1.0)
	θ = atan(agent.vel[2],agent.vel[1])
	cs, sn = siz.*(cos(θ), sin(θ))
	Polygon( Point2f.(map(x->[cs -sn; sn cs]*x, wejj)))
end

#-----------------------------------------------------------------------------------------
"""
	turn!( agent::AbstractAgent, angle)
	
Rotate agent's facing direction (vel) through the given angle.
"""
function turn!(agent::AbstractAgent, θ=pi)
	cs = cos(θ); sn = sin(θ)
	agent.vel = Tuple([cs -sn;sn cs]*collect(agent.vel))
end

"left!: Rotate agent's facing direction (vel) through a positive angle θ."
left!(agent::AbstractAgent, θ=pi/2) = turn!(agent, θ)

"right!: Rotate agent's facing direction (vel) through a negative angle θ."
right!(agent::AbstractAgent, θ=pi/2) = turn!(agent, -θ)

"wiggle!: Rotate agent's facing direction (vel) through a random angle of magnitude ≤ θ."
wiggle!(agent::AbstractAgent, θ=pi/4) = turn!(agent, (2rand()-1)*θ)

#-----------------------------------------------------------------------------------------
"""
	gradient( pos, F::AbstractArray, model::ABM)
	
Return the gradient vector for the given field F at this position in the model space.
"""
function gradient( pos, F::AbstractArray, model::ABM; h::Float64=1.0)
	nbhdfield = map([(-h,0),(h,0),(0,-h),(0,h)]) do step
		F[ get_spatial_index( normalize_position(pos.+step,model), F, model) ]
	end

	((nbhdfield[2]-nbhdfield[1]), (nbhdfield[4]-nbhdfield[3])) ./ 2h
end

#-----------------------------------------------------------------------------------------
"""
	valleys( extent::Tuple{Int,Int})

Return a 2D map of a multimodal valleys landscape with given extent.
"""
function valleys( extent::Tuple{Int,Int})
	xs = repeat( range(-3,3,extent[1]), 1, extent[2])
	ys = repeat( range(-3,3,extent[2])', extent[1], 1)

	map(
		(x,y) -> (1/3)*exp(-((x + 1)^2) - (y^2)) +
			10*(x/5 - (x^3) - (y^5)) * exp(-(x^2) - (y^2)) -
			(3*((1 - x)^2)) * exp(-(x^2) - ((y + 1)^2)),
		xs, ys
	)
end

#-----------------------------------------------------------------------------------------
"""
	dejong2( extent::Tuple{Int,Int})

Return a 2D map of a De Jong 2 landscape with given extent.
"""
function dejong2( extent::Tuple{Int,Int})
	xs = repeat( range(-10,10,extent[1]), 1, extent[2])
	ys = repeat( range(-10,10,extent[2])', extent[1], 1)

	map(
		(x,y) -> sin(2x) / (abs(x)+1) + sin(2y) / (abs(y)+1),
		xs, ys
	)
end

#-----------------------------------------------------------------------------------------
"""
	diffuse4( F::Matrix{R}, Λ=1.0) where R <: Real

Diffuse the field F via 4-neighbourhoods with periodic boundary conditions and the given
diffusion rate Λ.
"""
function diffuse4( F::Matrix{R}, Λ=1.0) where R <: Real
	F + 0.25Λ * (
		circshift(F,( 1,0)) + circshift(F,(0, 1)) +
		circshift(F,(-1,0)) + circshift(F,(0,-1)) - 4F
	)
end

#-----------------------------------------------------------------------------------------
"""
	diffuse4!( F::Matrix{R}, Λ=1.0) where R <: Real

Diffuse the field F in place via 4-neighbourhoods with periodic boundary conditions and the given
diffusion rate Λ.
"""
function diffuse4!( F::Matrix{R}, Λ=1.0) where R <: Real
	F[:] = diffuse4(F,Λ)[:]
	F
end

#-----------------------------------------------------------------------------------------
"""
	diffuse8( F::Matrix{R}, Λ=1.0) where R <: Real

Diffuse the field F via 8-neighbourhoods with periodic boundary conditions and the given
diffusion rate Λ.
"""
function diffuse8( F::Matrix{R}, Λ=1.0) where R <: Real
	sqrt2 = sqrt(2)						# Some preliminary calculations ...
	oneoversqrt2 = 0.5sqrt2
	sum_contributions = 4 + 2sqrt2

	F + (Λ/sum_contributions) * (
		circshift(F,( 1,0)) + circshift(F,( 0, 1)) +
		circshift(F,(-1,0)) + circshift(F,( 0,-1)) +
		oneoversqrt2*(
			circshift(F,( 1,1)) + circshift(F,( 1,-1)) +
			circshift(F,(-1,1)) + circshift(F,(-1,-1))
		) - (sum_contributions*F)
	)
end

#-----------------------------------------------------------------------------------------
"""
	diffuse8!( F::Matrix{R}, Λ=1.0) where R <: Real

Diffuse the field F in place via 8-neighbourhoods with periodic boundary conditions and the given
diffusion rate Λ.
"""
function diffuse8!( F::Matrix{R}, Λ=1.0) where R <: Real
	F[:] = diffuse8(F,Λ)[:]
	F
end

#-----------------------------------------------------------------------------------------
"""
	abmplayground( abm_constructor, kwargs...)

Extends abmexploration to replace the model of the ABMObservable abmobs when user clicks the
"reset model"-Button. This reproduces the behaviour of a NetLogo interface, in which resetting
leads to full reinitialisation of slider settings and agent population. (Nick Diercksen)
"""
function abmplayground( abm_constructor; kwargs...)
	playground, playobserver = Agents.abmexploration( abm_constructor(); kwargs...)

	# Locate Reset button from playground:
	reset_btn = nothing
	for element ∈ playground.content
		if element isa Button && occursin("reset", lowercase(element.label[]))
			reset_btn = element
			break
		end
	end
	@assert !isnothing(reset_btn) "Couldn't find the button 'Reset-model'!"

	# Retrieve initialisation keyword arguments from the abm_constructor function:
	init_kwargs = Base.kwarg_decl(@which abm_constructor())

	# Replace model with newly initialised one when reset button is clicked:
	on(reset_btn.clicks) do _
		props = abmproperties(playobserver.model.val)			# Retrieve current model properties
		kws = (; Dict([Pair(k, props[k]) for k in (keys(props) ∩ init_kwargs)])...)
		playobserver.model[] = abm_constructor(; kws...)		# Replace old model with new one
	end

	(playground,playobserver)
end

end # ... of module AgentTools
