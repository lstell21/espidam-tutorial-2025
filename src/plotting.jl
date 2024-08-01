"""
    plot_degree_distribution(graph_measures)

Plot the degree distribution of a graph.

# Arguments
- `graph_measures`: A structure containing graph measures.

# Returns
- `p`: A plot of the degree distribution.

# Example
```julia
plot_degree_distribution(graph_measures)
```
"""
function plot_degree_distribution(graph_measures::DataFrame)
    degree_distribution = graph_measures.degree_distribution[1]
    degree_distribution = DataFrame(degree=collect(keys(degree_distribution)), count=collect(values(degree_distribution)))
    p = plot(degree_distribution.degree, degree_distribution.count, seriestype=:bar, xlabel="Degree", ylabel="Count", legend=:none)
    return p
end

"""
    plot_epidemic_trajectories(mdf)

Plot the epidemic trajectories of susceptible, infected, and recovered individuals over time.

# Arguments
- `mdf`: A DataFrame containing the epidemic data with columns `:susceptible_count`, `:infected_count`, and `:recovered_count`.

# Example
```julia
plot_epidemic_trajectories(mdf)
```
"""
function plot_epidemic_trajectories(mdf::DataFrame)
    # Extract the data from mdata
    susceptible = mdf[!, :susceptible_count]
    infected = mdf[!, :infected_count]
    recovered = mdf[!, :recovered_count]

    # Create a time vector
    time = 1:length(susceptible)

    # Plot the trajectories
    plot(time, susceptible, label="Susceptible", legend=:topright, xlabel="Time", ylabel="Count", linewidth=2)
    plot!(time, infected, label="Infected", legend=:topright, linewidth=2)
    plot!(time, recovered, label="Recovered", legend=:topright, linewidth=2)
end

"""
    plot_single_run(network_type::Symbol, mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1, patient_zero::Symbol=:random, high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, n_steps::Int=100)

Plot a single run of an epidemic simulation.

Arguments:
- `network_type::Symbol`: The type of network to use for the simulation. Possible values are `:random`, `:smallworld`, `:preferentialattachment`, `:configuration`, or `:proportionatemixing`.
- `mean_degree::Int`: The mean degree of the network. Default is 4.
- `n_nodes::Int`: The number of nodes in the network. For :configuration, the number of nodes is fixed to 1000. Default is 1000. 
- `dispersion::Float64`: The dispersion parameter for the network. Default is 0.1.
- `patient_zero::Symbol`: The type of patient zero to use for the simulation. Can be `:random`(a random agent), `:maxdegree` (the agent with highest degree_centrality), `:maxbetweenness`, and `:maxeigenvector`. Default is `:random`.
- `high_risk::Symbol`: How high-risk individuals are distributed in the graph. Can be `:random` (randomly distributed), `:maxdegree` (based on degree centrality), `:maxbetweenness` (based on betweenness centrality), and `:maxeigenvector` (based on eigenvector centrality). Default is `:random`.
- `fraction_high_risk::Float64`: The fraction of high-risk individuals in the population. Default is 0.1.
- `trans_prob::Float64`: The transmission probability. Default is 0.1.
- `n_steps::Int`: The number of simulation steps to run. Default is 100.

Returns:
- Nothing.

Example:
```julia
plot_single_run(:random, :random, :random, 4, 100)
```
"""
function plot_single_run(; network_type::Symbol,  mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1, patient_zero::Symbol=:random, high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, n_steps::Int=100)
    model = initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob)
    _, mdf = run!(model, n_steps; adata, mdata)
    graph_measures = analyze_graph(model.graph)
    plotdynamics=plot_epidemic_trajectories(mdf)
    savefig(plotdynamics, "figures/plotdynamics_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).pdf")
    plotdegdist = plot_degree_distribution(graph_measures)
    display(plotdegdist)
    savefig(plotdegdist,"figures/plotdegdist_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).pdf")
end


"""
    run_and_plot(plot_type::Symbol, network_type::Symbol, mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1, patient_zero::Symbol=:random,  high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, n_steps::Int=100)

Run simulations, save the results to file and plot the results based on the specified parameters.

# Arguments
- `plot_type::Symbol`: The type of plot to generate. Possible values are `:sfr` (susceptible fraction remaining), `:max_infected` (maximum number of infected), or `:duration` (duration of epidemic).
- `network_type::Symbol`: The type of network to use for the simulation. Possible values are `:random`, `:smallworld`, `:preferentialattachment`, `:configuration`, or `:proportionatemixing`.
- `mean_degree::Int`: The mean degree of the network. Default is 4.
- `n_nodes::Int`: The number of nodes in the network. For :configuration, the number of nodes is fixed to 1000. Default is 1000.
- `dispersion::Float64`: The dispersion parameter for the network. Default is 0.1.
- `patient_zero::Symbol`: The type of patient zero to use for the simulation. Can be `:random`(a random agent), `:maxdegree` (the agent with highest degree_centrality), `:maxbetweenness`, and `:maxeigenvector`. Default is `:random`.
- `high_risk::Symbol`: How high-risk individuals are distributed in the graph. Can be `:random` (randomly distributed), `:maxdegree` (based on degree centrality), `:maxbetweenness` (based on betweenness centrality), and `:maxeigenvector` (based on eigenvector centrality). Default is `:random`.
- `fraction_high_risk::Float64`: The fraction of high-risk individuals in the population. Default is `0.1`.
- `trans_prob::Float64`: The transmission probability. Default is 0.1.
- `n_steps::Int`: The number of simulation steps to run. Default is 100.

# Example
```julia
run_and_plot(:sfr, :proportionatemixing, :random, 4)
```
"""
function run_and_plot(; plot_type::Symbol, network_type::Symbol,  mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1, patient_zero::Symbol=:random, high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, n_steps::Int=100)
    model = initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob)
    graph_measures = analyze_graph(model.graph)
    CSV.write("data/graph_keyfig_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).csv", graph_measures)
    multiple_runs = run_simulations(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob, n_steps)
    grouped_data = groupby(multiple_runs, [:seed])
    final_results = combine(grouped_data, :infected_count => argmin => :first_to_last_infected, :infected_count => maximum => :max_infected)
    last_rows = combine(grouped_data, names(multiple_runs) .=> last)
    final_results[!, :susceptible_fraction_remaining] = last_rows.susceptible_count_last ./ (last_rows.susceptible_count_last + last_rows.infected_count_last + last_rows.recovered_count_last)
    # Write the results to a CSV file
    CSV.write("data/simulation_results_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).csv", multiple_runs)
    CSV.write("output/final_results_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).csv", final_results)
    if plot_type == :sfr
        plot(final_results.susceptible_fraction_remaining, seriestype=:bar, xlabel="Seeds", ylabel="Susceptible Fraction Remaining", legend=:none)
    elseif plot_type == :max_infected
        plot(final_results.max_infected, seriestype=:bar, xlabel="Seeds", ylabel="Maximum Number of Infected", legend=:none)
    elseif plot_type == :duration
        plot(final_results.first_to_last_infected, seriestype=:bar, xlabel="Seeds", ylabel="Duration of Epidemic", legend=:none)
    end
end