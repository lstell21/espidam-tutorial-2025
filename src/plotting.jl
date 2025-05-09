"""
    plot_degree_distribution(degree_distribution; network_type::Symbol)

Plot the degree distribution of a graph.

# Arguments
- `degree_distribution`: A dictionary mapping degrees to counts or DataFrame containing degree distribution
- `network_type`: Symbol indicating the type of network (for title)

# Returns
- `p`: A plot of the degree distribution.

# Example
```julia
plot_degree_distribution(graph_analysis["degree_distribution"], :random)
```
"""
function plot_degree_distribution(degree_distribution; network_type::Symbol)
    # Handle the case when degree_distribution is a dictionary (from analyze_graph)
    if isa(degree_distribution, Dict)
        degrees = collect(keys(degree_distribution))
        counts = collect(values(degree_distribution))
    # Handle the case when degree_distribution is a DataFrame with degree_distribution column
    elseif isa(degree_distribution, DataFrame) && haskey(degree_distribution, :degree_distribution)
        dd = degree_distribution.degree_distribution[1]
        degrees = collect(keys(dd))
        counts = collect(values(dd))
    # Handle the case when degree_distribution is already unpacked DataFrame
    elseif isa(degree_distribution, DataFrame) && haskey(degree_distribution, :degree) && haskey(degree_distribution, :count)
        degrees = degree_distribution.degree
        counts = degree_distribution.count
    else
        error("Unsupported degree_distribution format")
    end
    
    p = plot(degrees, counts, seriestype=:bar, title="Degree Distribution ($(String(network_type)))", 
             xlabel="Degree", ylabel="Count", legend=:none)
    return p
end

"""
    plot_epidemic_trajectories(mdf, network_type::Symbol)

Plot the epidemic trajectories of susceptible, infected, and recovered individuals over time.

# Arguments
- `mdf`: A DataFrame containing the epidemic data with columns `:susceptible_count`, `:infected_count`, and `:recovered_count`.
- `network_type`: Symbol indicating the type of network (for title).

# Returns
- `p`: A plot of the epidemic trajectories.

# Example
```julia
plot_epidemic_trajectories(mdf, :random)
```
"""
function plot_epidemic_trajectories(mdf::DataFrame, network_type::Symbol)
    # Extract the data from mdata
    susceptible = mdf[!, :susceptible_count]
    infected = mdf[!, :infected_count]
    recovered = mdf[!, :recovered_count]

    # Create a time vector
    time = 1:length(susceptible)

    # Plot the trajectories with network type in the title
    p = plot(time, susceptible, label="Susceptible", legend=:topright, 
             xlabel="Time", ylabel="Count", linewidth=2,
             title="Epidemic Dynamics ($(String(network_type)))")
    plot!(p, time, infected, label="Infected", linewidth=2)
    plot!(p, time, recovered, label="Recovered", linewidth=2)

    return p
end

"""
    plot_single_run(; network_type::Symbol, mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1, patient_zero::Symbol=:random, high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, n_steps::Int=100)

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
- `plotdynamics`: A plot of the epidemic trajectories (susceptible, infected, recovered) over time.
- `plotdegdist`: A plot of the degree distribution of the graph.
- `combined_plot`: A combined plot with both the epidemic trajectories and degree distribution side-by-side.

Example:
```julia
dynamics_plot, degdist_plot, combined_plot = plot_single_run(network_type=:random)
```
"""
function plot_single_run(; network_type::Symbol,  mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1, patient_zero::Symbol=:random, high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, n_steps::Int=100)
    model = initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob)
    
    # Define adata and mdata locally to avoid relying on global variables
    adata = [:status]
    mdata = [:susceptible_count, :infected_count, :recovered_count]
    
    _, mdf = run!(model, n_steps; adata, mdata)
    plotdynamics = plot_epidemic_trajectories(mdf, model.network_type)
    savefig(plotdynamics, "figures/plotdynamics_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).pdf")
    
    # Get graph analysis results
    graph_analysis = analyze_graph(model.graph)
    
    # Plot degree distribution - pass the dictionary directly
    plotdegdist = plot_degree_distribution(graph_analysis["degree_distribution"]; network_type=model.network_type)
    savefig(plotdegdist,"figures/plotdegdist_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).pdf")
    
    # Create a combined plot for side-by-side visualization with network type in titles
    combined_plot = plot(plotdynamics, plotdegdist, layout=(1,2), size=(1000, 400),
                         title=["Epidemic Dynamics ($(String(model.network_type)))" "Degree Distribution ($(String(model.network_type)))"])
    savefig(combined_plot,"figures/combined_plot_$(model.network_type)_mdeg_$(model.mean_degree).pdf")
    
    # Return all three plots
    return plotdynamics, plotdegdist, combined_plot
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

# Returns
- A plot based on the specified plot_type

# Example
```julia
result_plot = run_and_plot(plot_type=:sfr, network_type=:proportionatemixing, patient_zero=:random)
```
"""
function run_and_plot(; plot_type::Symbol, network_type::Symbol,  mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1, patient_zero::Symbol=:random, high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, n_steps::Int=100)
    model = initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob)
    multiple_runs = run_simulations(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob, n_steps)
    grouped_data = groupby(multiple_runs, [:seed])
    final_results = combine(grouped_data, :infected_count => argmin => :first_to_last_infected, :infected_count => maximum => :max_infected)
    last_rows = combine(grouped_data, names(multiple_runs) .=> last)
    final_results[!, :susceptible_fraction_remaining] = last_rows.susceptible_count_last ./ (last_rows.susceptible_count_last + last_rows.infected_count_last + last_rows.recovered_count_last)
    # Write the results to a CSV file
    CSV.write("data/simulation_results_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).csv", multiple_runs)
    CSV.write("output/final_results_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).csv", final_results)
    
    # Create descriptive title
    network_desc = "$(titlecase(String(network_type))) Network (mean degree: $(mean_degree))"
    patient_desc = patient_zero == :random ? "Random Patient Zero" : "Patient Zero: $(titlecase(String(patient_zero)))"
    plot_title = ""
    
    plot_title = "$(plot_type == :sfr ? "Susceptible Fraction Remaining" : plot_type == :max_infected ? "Maximum Number of Infected" : "Epidemic Duration")\n$(network_desc)"
    ylabel = plot_type == :sfr ? "Susceptible Fraction Remaining" : plot_type == :max_infected ? "Maximum Number of Infected" : "Duration of Epidemic (steps)"
    data = plot_type == :sfr ? final_results.susceptible_fraction_remaining : plot_type == :max_infected ? final_results.max_infected : final_results.first_to_last_infected

    result_plot = plot(data, 
                       seriestype=:bar, 
                       xlabel="Simulation Run", 
                       ylabel=ylabel, 
                       legend=:none,
                       title=plot_title,
                       guidefontsize=10,
                       titlefontsize=12,
                       size=(700, 500),
                       margin=10mm)
    
    # Save the plot with a descriptive filename
    savefig(result_plot, "figures/$(plot_type)_plot_$(model.network_type)_mdeg_$(model.mean_degree).pdf")
    
    # Return the plot for display in notebook
    return result_plot
end