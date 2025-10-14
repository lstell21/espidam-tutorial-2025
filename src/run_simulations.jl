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
- `low_risk_factor::Float64`: Factor to multiply the transmission probability for low risk agents. Must be between 0 and 1. Default is 1.0.

# Returns
- `mdf::DataFrame`: A DataFrame containing the simulation results.

# Example
```julia
mdf = run_simulations(network_type=:random, mean_degree=4, patient_zero=:random, high_risk=:random, fraction_high_risk=1.0)
```
"""
function run_simulations(; network_type::Symbol, mean_degree::Int, n_nodes::Int=1000, 
                        dispersion::Float64=0.1, patient_zero::Symbol=:random, 
                        high_risk::Symbol=:random, fraction_high_risk::Float64=1.0, low_risk_factor::Float64=1.0,
                        trans_prob::Float64=0.1, n_steps::Int=100, r̂=nothing, p̂=nothing)
    # Validate low_risk_factor
    if !(0 <= low_risk_factor <= 1)
        error("low_risk_factor must be between 0 and 1, got $low_risk_factor")
    end
    
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
        :low_risk_factor => low_risk_factor,
        :days_to_recovered => 14
    )
    
    # Add r̂ and p̂ to parameters if provided
    if r̂ !== nothing
        parameters[:r̂] = r̂
    end
    if p̂ !== nothing
        parameters[:p̂] = p̂
    end

    # Data to collect
    adata = [:status]
    mdata = [:susceptible_count, :infected_count, :hospitalized_count, :recovered_count]
    
    # Run the simulation for all combinations of parameters
    _, mdf = paramscan(
        parameters,
        initialize;
        mdata=mdata,
        n=n_steps,
        showprogress=false
    );
    
    # Calculate summary statistics
    summary_stats = DataFrame()
    for step in 0:n_steps
        step_data = filter(row -> row.step == step, all_mdata)
        if nrow(step_data) > 0
            summary_row = DataFrame(
                step = step,
                susceptible_count_mean = mean(step_data.susceptible_count),
                susceptible_count_std = std(step_data.susceptible_count),
                infected_count_mean = mean(step_data.infected_count),
                infected_count_std = std(step_data.infected_count),
                hospitalized_count_mean = mean(step_data.hospitalized_count),
                hospitalized_count_std = std(step_data.hospitalized_count),
                recovered_count_mean = mean(step_data.recovered_count),
                recovered_count_std = std(step_data.recovered_count)
            )
            
            # Add confidence intervals (95%)
            n_sims = nrow(step_data)
            summary_row.infected_count_ci_lower = summary_row.infected_count_mean[1] - 1.96 * summary_row.infected_count_std[1] / sqrt(n_sims)
            summary_row.infected_count_ci_upper = summary_row.infected_count_mean[1] + 1.96 * summary_row.infected_count_std[1] / sqrt(n_sims)
            summary_row.hospitalized_count_ci_lower = summary_row.hospitalized_count_mean[1] - 1.96 * summary_row.hospitalized_count_std[1] / sqrt(n_sims)
            summary_row.hospitalized_count_ci_upper = summary_row.hospitalized_count_mean[1] + 1.96 * summary_row.hospitalized_count_std[1] / sqrt(n_sims)
            
            append!(summary_stats, summary_row)
        end
    end
    
    return summary_stats
end