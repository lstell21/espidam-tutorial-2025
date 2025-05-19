"""
    run_simulations(network_type::Symbol, mean_degree::Int, patient_zero::Symbol, high_risk::Symbol, fraction_high_risk::Float64)

Run simulations for an epidemiological model with different combinations of parameters.

# Arguments
- `network_type::Symbol`: The type of network to use for the simulation.
- `patient_zero::Symbol`: The type of patient zero to use for the simulation. Default is `:random`.
- `mean_degree::Int`: The mean degree of the network. Default is 4.
- `n_nodes::Int`: The number of nodes in the network. Default is 1000.
- `dispersion::Float64`: The dispersion of the network. Default is 0.1.
- `high_risk::Symbol`: The type of high-risk individuals to consider. Default is `:random`.
- `fraction_high_risk::Float64`: The fraction of high-risk individuals in the population. Default is 1.0.
- `trans_prob::Float64`: The transmission probability. Default is 0.1.
- `n_steps::Int`: The number of simulation steps to run. Default is 100.
- `r̂`: The r parameter for negative binomial distribution, used only when `network_type` is `:proportionatemixing`. Default is nothing.
- `p̂`: The p parameter for negative binomial distribution, used only when `network_type` is `:proportionatemixing`. Default is nothing.

# Returns
- `mdf::DataFrame`: A DataFrame containing the simulation results.

# Example
```julia
mdf = run_simulations(network_type=:random, mean_degree=4, patient_zero=:random, high_risk=:random, fraction_high_risk=0.1)
```
"""
function run_simulations(; network_type::Symbol, mean_degree::Int, n_nodes::Int=1000, 
                        dispersion::Float64=0.1, patient_zero::Symbol=:random, 
                        high_risk::Symbol=:random, fraction_high_risk::Float64=1.0, 
                        trans_prob::Float64=0.1, n_steps::Int=100, r̂=nothing, p̂=nothing)
    # Define parameters
    parameters = Dict(
        :seed => rand(UInt16, 100),
        :network_type => network_type,
        :mean_degree => mean_degree,
        :n_nodes => n_nodes,
        :dispersion => dispersion,
        :patient_zero => patient_zero,
        :high_risk => high_risk,
        :trans_prob => trans_prob,
        :fraction_high_risk => fraction_high_risk,
        :days_to_recovered => 14
    )
    
    # Add r̂ and p̂ to parameters if provided
    if r̂ !== nothing
        parameters[:r̂] = r̂
    end
    if p̂ !== nothing
        parameters[:p̂] = p̂
    end

    # Run the simulation for all combinations of parameters
    _, mdf = paramscan(
        parameters,
        initialize;
        mdata=[:susceptible_count, :infected_count, :recovered_count],
        n=n_steps,
        showprogress=false
    );
    return mdf
end