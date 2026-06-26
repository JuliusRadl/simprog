#========================================================================================#
"""
	DAH

Implement a demonstration of Newman and Müller's Differential Adhesion Hypothesis (DAH).

Author: Niall Palfreyman, March 2025.
"""
module DAH

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	Cell

A Cell is an Agent possessing the property `cadherin`, which specifies the stickiness of the cell
to other cells in its neighbourhood.
"""
@agent struct Cell(ContinuousAgent{2,Float64})
	cadherin::Float64
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	dah( kwargs)

Initialise the DAH model.
"""
function dah(;
	adhesion_range=4.0,		# Max separation within which two Cells can adhere
	adhesion=0.5,			# Base adhesion level between interacting Cells
	radius=1.2,				# Incompressibility radius of Cells
	repulsion=0.5,			# Repulsive force between incompressible Cells
	thermal=0.01,			# Thermal force on Cells
	gravity=0.0,			# Gravitational force on Cells (very small, e.g. 0.01)
	adherins=3				# Number of distinct adherin classes
)
	properties = Dict(
		:adhesion_range => adhesion_range,
		:adhesion => adhesion,
		:radius => radius,
		:repulsion => repulsion,
		:thermal => thermal,
		:gravity => gravity,
		:adherins => adherins,
	)

	extent = (20,20)			# Extent of world
	agents_per_square = 3		# Number of agents per unit square in the world

	dah = StandardABM( Cell, ContinuousSpace(extent;spacing=1.0);
		agent_step!,
		properties,
	)

	for _ in 1:agents_per_square*prod(extent)
		add_agent!( dah, (0,0), rand(1:adherins))
	end

	return dah
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!(cell,dah)

Each Cell in the DAH model moves under the influence of three forces: adhere to a random nearby
cell, repel a nearby cell if it is too close, meander under the influence of thermal motion,
and possibly fall towards the world's origin.
"""
function agent_step!( cell::Cell, dah::ABM)
	adhere!( cell, dah)
	repel!( cell, dah)
	meander!( cell, dah)
	fall!( cell, dah)
end

#-----------------------------------------------------------------------------------------
"""
	adhere!(cell,dah)

Adhere towards a neighbouring cell within a radius of adhesion_range
"""
function  adhere!( cell::Cell, dah::ABM)
	nbr = random_nearby_agent( cell, dah, dah.adhesion_range)
	if nbr !== nothing
		cell.vel = get_direction( cell.pos, nbr.pos, dah)
		move_agent!( cell, dah,
			dah.adhesion * stickiness(cell.cadherin,nbr.cadherin,dah.adherins)
		)
	end
end

#-----------------------------------------------------------------------------------------
"""
	stickiness(cad1,cad2,nadherins)

Calculate the cadherin stickiness between two cadherin classes, assuming the given total number
of possible adherin classes. This formula is provisional and VERY suspicious!
"""
function  stickiness( cad1, cad2, nadherins)
	1 - (2abs(cad1 - cad2) + cad1 + cad2) / 2nadherins
end

#-----------------------------------------------------------------------------------------
"""
	repel!(cell,dah)

Repel myself away from a nearby cell that has approached closer than my radius.
"""
function  repel!( cell::Cell, dah::ABM)
	nbr = random_nearby_agent( cell, dah, dah.radius)
	if nbr !== nothing
		cell.vel = get_direction( cell.pos, nbr.pos, dah)
		move_agent!( cell, dah, -dah.repulsion)
	end       
end

#-----------------------------------------------------------------------------------------
"""
	meander!(cell,dah)

Thermal motion: Turn in a random direction and wander a short way in that direction.
"""
function  meander!( cell::Cell, dah::ABM)
	turn!(cell,2π*rand())
	move_agent!( cell, dah, dah.thermal)
end

#-----------------------------------------------------------------------------------------
"""
	fall!(cell,dah)

Gravitational motion: Gravitate towards the world's centre: abmspace(dah).extent./2).
??? Make this stronger!
"""
function fall!( cell::Cell, dah::ABM)
	if dah.gravity > 0.0
#		cell.vel = get_direction( cell.pos, abmspace(dah).extent./2, dah)
		cell.vel = (0, abmspace(dah).extent[2]/2-cell.pos[2])
		move_agent!( cell, dah, dah.gravity)
	end
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Set up a playground for differential adhesion.
"""
function demo()
	params = Dict(
		:adhesion_range => 1:7,
		:adhesion => 0:0.1:1,
		:radius => 0.5:0.1:2,
		:repulsion => 0:0.1:1,
		:thermal => 0:0.001:0.1,
		:gravity => 0:0.01:1,
		:adherins => 1:length(spectrum),
	)
	plotkwargs = (
		agent_color = (ag->spectrum[Int(ag.cadherin)]),
		agent_size = 25,
	)

	playground, = abmplayground( dah; params, plotkwargs...)

	playground
end

end