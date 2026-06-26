#========================================================================================#
"""
	Turing

This simulation demonstrates the emergence of complex activity in reaction-diffusion
(i.e., Turing) systems.

Author: Niall Palfreyman, March 2025
"""
module Turing

include( "../../Development/Generative/AgentTools.jl")

using Agents, GLMakie, Random, .AgentTools

#-----------------------------------------------------------------------------------------
# Module types:
#-----------------------------------------------------------------------------------------
"""
	SlimeMould

A SlimeMould can in principle move, and secretes the chemical A(ctivator).
"""
@agent struct SlimeMould(ContinuousAgent{2,Float64})
end

#-----------------------------------------------------------------------------------------
# Module methods:
#-----------------------------------------------------------------------------------------
"""
	turing( kwargs)

Initialise the Turing reaction-diffusion system.
"""
function turing(;
	a_secr_rate = 0.0,
	a_diff_rate = 1e-2,
	i_diff_rate = 0.4,
	a_evap_rate = 1e-2,
	i_evap_rate = 3e-2,
	pPop = 0.1,
)
	extent = (50,50)
	properties = Dict(
		:activator => 1e-4rand(extent...),
		:inhibitor => 1e-4rand(extent...),
		:a_secr_rate => a_secr_rate,
		:a_diff_rate => a_diff_rate,
		:i_diff_rate => i_diff_rate,
		:a_evap_rate => a_evap_rate,
		:i_evap_rate => i_evap_rate,
		:pPop => pPop,
		:dt => 0.5,
	)

	tur = StandardABM( SlimeMould, ContinuousSpace(extent,spacing=1.0);
		agent_step!, model_step!, properties
	)

	for _ in 1:pPop*prod(extent)
		add_agent!( SlimeMould, tur, (0,0))
	end

	tur
end

#-----------------------------------------------------------------------------------------
"""
	agent_step!( slimy, tur)

Express activator and inhibitor proteins.
"""
function agent_step!( slimy::SlimeMould, tur::ABM)
	idx_here = get_spatial_index( slimy.pos, tur.activator, tur)
	if tur.a_secr_rate > 0.0
		tur.activator[idx_here] += tur.a_secr_rate*tur.dt
	end
end

#-----------------------------------------------------------------------------------------
"""
	model_step!(tur)

Allow both A and I to react, evaporate and diffuse.
"""
function model_step!( tur::ABM)
	# Chemical reaction:
	(tur.activator,tur.inhibitor) = (
		tur.activator + tur.dt * tur.activator.^2 ./ tur.inhibitor,
		tur.inhibitor + tur.dt * tur.activator.^2
	)
	
	# Evaporation:
	tur.activator *= (1 - tur.dt*tur.a_evap_rate)
	tur.inhibitor *= (1 - tur.dt*tur.i_evap_rate)

	# Diffusion:
	diffuse4!( tur.activator, tur.dt*tur.a_diff_rate)
	diffuse4!( tur.inhibitor, tur.dt*tur.i_diff_rate)
end

#-----------------------------------------------------------------------------------------
"""
	demo()

Create playground with two wave sources.
"""
function demo()
	params = Dict(
		:pPop => 0:0.01:1,
		:a_secr_rate => 0:0.001:0.01,
		:a_diff_rate => 0:0.001:0.025,
		:i_diff_rate => 0:0.01:1,
		:a_evap_rate => 0:0.001:0.1,
		:i_evap_rate => 0:0.001:0.1,
	)
	plotkwargs = (
		agent_marker = :circle,
		agent_color = :red,
		agent_size = 15,
		heatarray = :activator,
		add_colorbar = false,
	)

	playground, = abmplayground( turing; params, plotkwargs...)

	playground
end

end