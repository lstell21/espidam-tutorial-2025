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
    run_and_plot_comparison(; network_types::Vector{Symbol}, mean_degree::Int=4, n_nodes::Int=1000, 
                          dispersion::Float64=0.1, patient_zero::Symbol=:random, high_risk::Symbol=:random, 
                          fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, n_steps::Int=100)

Run simulations for multiple network types, save the results to file and generate comparison plots 
of epidemic metrics across different network types.

# Arguments
- `network_types::Vector{Symbol}`: Vector of network types to compare. Possible values include `:random`, 
  `:smallworld`, `:preferentialattachment`, `:configuration`, or `:proportionatemixing`.
- `mean_degree::Int`: The mean degree of the network. Default is 4.
- `n_nodes::Int`: The number of nodes in the network. Default is 1000.
- `dispersion::Float64`: The dispersion parameter for the network. Default is 0.1.
- `patient_zero::Symbol`: The type of patient zero to use for the simulation. Default is `:random`.
- `high_risk::Symbol`: How high-risk individuals are distributed. Default is `:random`.
- `fraction_high_risk::Float64`: The fraction of high-risk individuals in the population. Default is `0.1`.
- `trans_prob::Float64`: The transmission probability. Default is 0.1.
- `n_steps::Int`: The number of simulation steps to run. Default is 100.

# Returns
- A tuple containing three comparison plots: (duration_comparison, max_infected_comparison, sfr_comparison)

# Example
```julia
duration_comparison, max_infected_comparison, sfr_comparison = run_and_plot_comparison(
    network_types=[:random, :smallworld, :preferentialattachment],
    mean_degree=4
)
```
"""
function run_and_plot_comparison(; network_types::Vector{Symbol}, mean_degree::Int=4, n_nodes::Int=1000, 
                               dispersion::Float64=0.1, patient_zero::Symbol=:random, high_risk::Symbol=:random, 
                               fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, n_steps::Int=100)
    # Store results for each network type
    all_results = Dict()
    
    # Run simulations for each network type
    for network_type in network_types
        println("Running simulations for $(network_type) network...")
        model = initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob)
        multiple_runs = run_simulations(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob, n_steps);
        grouped_data = groupby(multiple_runs, [:seed])
        final_results = combine(grouped_data, :infected_count => argmin => :first_to_last_infected, :infected_count => maximum => :max_infected)
        last_rows = combine(grouped_data, names(multiple_runs) .=> last)
        final_results[!, :susceptible_fraction_remaining] = last_rows.susceptible_count_last ./ (last_rows.susceptible_count_last + last_rows.infected_count_last + last_rows.recovered_count_last)
        
        # Save results for each network type
        CSV.write("data/simulation_results_$(network_type)_mdeg_$(mean_degree)_nn_$(n_nodes)_disp_$(dispersion)_pat0_$(patient_zero)_hirisk_$(high_risk)_hr_frac_$(fraction_high_risk)_trans_$(trans_prob).csv", multiple_runs)
        CSV.write("output/final_results_$(network_type)_mdeg_$(mean_degree)_nn_$(n_nodes)_disp_$(dispersion)_pat0_$(patient_zero)_hirisk_$(high_risk)_hr_frac_$(fraction_high_risk)_trans_$(trans_prob).csv", final_results)
        
        # Calculate summary statistics for each network type
        all_results[network_type] = (
            duration_mean = mean(final_results.first_to_last_infected),
            duration_std = std(final_results.first_to_last_infected),
            max_infected_mean = mean(final_results.max_infected),
            max_infected_std = std(final_results.max_infected),
            sfr_mean = mean(final_results.susceptible_fraction_remaining),
            sfr_std = std(final_results.susceptible_fraction_remaining)
        )
    end
    
    # Prepare data for comparison plots
    network_labels = [titlecase(String(nt)) for nt in network_types]
    duration_means = [all_results[nt].duration_mean for nt in network_types]
    duration_stds = [all_results[nt].duration_std for nt in network_types]
    max_infected_means = [all_results[nt].max_infected_mean for nt in network_types]
    max_infected_stds = [all_results[nt].max_infected_std for nt in network_types]
    sfr_means = [all_results[nt].sfr_mean for nt in network_types]
    sfr_stds = [all_results[nt].sfr_std for nt in network_types]
    
    # Generate comparison plots with optimized layout
    
    # 1. Duration comparison
    duration_comparison = bar(
        network_labels, 
        duration_means,
        yerror=duration_stds,
        title="Epidemic Duration Comparison\n(mean degree: $(mean_degree))",
        xlabel="Network Type",
        ylabel="Duration (steps)",
        legend=false,
        size=(600, 300),         # Reduced size
        margin=5mm,              # Reduced margins
        guidefontsize=9,         # Smaller font size
        titlefontsize=10,        # Smaller title font
        xtickfontsize=8,         # Smaller tick labels
        bar_width=0.4,           # Thinner bars (reduced from 0.6)
        linewidth=1,             # Thinner error bars
        left_margin=8mm          # Add left margin for aligned y-axis label
    )
    savefig(duration_comparison, "figures/duration_comparison_mdeg_$(mean_degree).pdf")
    
    # 2. Maximum infected comparison
    max_infected_comparison = bar(
        network_labels,
        max_infected_means,
        yerror=max_infected_stds,
        title="Maximum Infected Comparison\n(mean degree: $(mean_degree))",
        xlabel="Network Type",
        ylabel="Maximum Number of Infected",
        legend=false,
        size=(600, 300),         # Reduced size
        margin=5mm,              # Reduced margins
        guidefontsize=9,         # Smaller font size
        titlefontsize=10,        # Smaller title font
        xtickfontsize=8,         # Smaller tick labels
        bar_width=0.4,           # Thinner bars (reduced from 0.6)
        linewidth=1,             # Thinner error bars
        left_margin=8mm          # Add left margin for aligned y-axis label
    )
    savefig(max_infected_comparison, "figures/max_infected_comparison_mdeg_$(mean_degree).pdf")
    
    # 3. Susceptible fraction remaining comparison
    sfr_comparison = bar(
        network_labels,
        sfr_means,
        yerror=sfr_stds,
        title="Susceptible Fraction Remaining Comparison\n(mean degree: $(mean_degree))",
        xlabel="Network Type",
        ylabel="Susceptible Fraction Remaining",
        legend=false,
        size=(600, 300),         # Reduced size
        margin=5mm,              # Reduced margins
        guidefontsize=9,         # Smaller font size
        titlefontsize=10,        # Smaller title font
        xtickfontsize=8,         # Smaller tick labels
        bar_width=0.4,           # Thinner bars (reduced from 0.6)
        linewidth=1,             # Thinner error bars
        ylims=(0, maximum(sfr_means + sfr_stds) * 1.1), # Adjust y-axis to better show small values
        left_margin=8mm          # Add left margin for aligned y-axis label
    )
    savefig(sfr_comparison, "figures/sfr_comparison_mdeg_$(mean_degree).pdf")
    
    # Create a more compact combined comparison plot with aligned y-axis labels
    combined_comparison = plot(
        duration_comparison, max_infected_comparison, sfr_comparison,
        layout=(3,1),
        size=(700, 700),         # More compact size
        margin=3mm,              # Smaller margins around the overall plot
        title=["Epidemic Duration" "Maximum Infected" "Susceptible Fraction Remaining"],
        titlefontsize=10,
        plot_title="Network Comparison (mean degree: $(mean_degree))", 
        plot_titlefontsize=12,
        left_margin=8mm          # Consistent left margin for alignment
    )
    savefig(combined_comparison, "figures/combined_comparison_mdeg_$(mean_degree).pdf")
    
    # Return only the combined comparison plot
    return combined_comparison
end


"""
    run_and_plot(; network_type::Symbol, mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1, patient_zero::Symbol=:random, high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, n_steps::Int=100)

Run simulations, save the results to file and generate three plots: duration of epidemic, maximum number of infected, and susceptible fraction remaining.

# Arguments
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
- A tuple containing all three plots: (duration_plot, max_infected_plot, sfr_plot)

# Example
```julia
duration_plot, max_infected_plot, sfr_plot = run_and_plot(network_type=:proportionatemixing, patient_zero=:random)
```
"""
function run_and_plot(; network_type::Symbol, mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1, patient_zero::Symbol=:random, high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, n_steps::Int=100)
    model = initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob)
    multiple_runs = run_simulations(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob, n_steps)
    grouped_data = groupby(multiple_runs, [:seed])
    final_results = combine(grouped_data, :infected_count => argmin => :first_to_last_infected, :infected_count => maximum => :max_infected)
    last_rows = combine(grouped_data, names(multiple_runs) .=> last)
    final_results[!, :susceptible_fraction_remaining] = last_rows.susceptible_count_last ./ (last_rows.susceptible_count_last + last_rows.infected_count_last + last_rows.recovered_count_last)
    
    # Write the results to CSV files
    CSV.write("data/simulation_results_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).csv", multiple_runs)
    CSV.write("output/final_results_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).csv", final_results)
    
    # Create descriptive base title that includes network type
    network_desc = "$(titlecase(String(network_type))) Network (mean degree: $(mean_degree))"
    
    # Generate the three plots
    
    # 1. Duration plot
    duration_plot = plot(final_results.first_to_last_infected, 
                     seriestype=:bar, 
                     xlabel="Simulation Run", 
                     ylabel="Duration of Epidemic (steps)",
                     legend=:none,
                     title="Epidemic Duration\n$(network_desc)",
                     guidefontsize=10,
                     titlefontsize=12,
                     size=(700, 500),
                     margin=10mm)
    savefig(duration_plot, "figures/duration_plot_$(model.network_type)_mdeg_$(model.mean_degree).pdf")
    
    # 2. Maximum infected plot
    max_infected_plot = plot(final_results.max_infected, 
                         seriestype=:bar, 
                         xlabel="Simulation Run", 
                         ylabel="Maximum Number of Infected",
                         legend=:none,
                         title="Maximum Number of Infected\n$(network_desc)",
                         guidefontsize=10,
                         titlefontsize=12,
                         size=(700, 500),
                         margin=10mm)
    savefig(max_infected_plot, "figures/max_infected_plot_$(model.network_type)_mdeg_$(model.mean_degree).pdf")
    
    # 3. Susceptible fraction remaining plot
    sfr_plot = plot(final_results.susceptible_fraction_remaining, 
                seriestype=:bar, 
                xlabel="Simulation Run", 
                ylabel="Susceptible Fraction Remaining",
                legend=:none,
                title="Susceptible Fraction Remaining\n$(network_desc)",
                guidefontsize=10,
                titlefontsize=12,
                size=(700, 500),
                margin=10mm)
    savefig(sfr_plot, "figures/sfr_plot_$(model.network_type)_mdeg_$(model.mean_degree).pdf")
    
    # Create a combined plot with all three metrics
    combined_plot = plot(duration_plot, max_infected_plot, sfr_plot, 
                    layout=(3,1), 
                    size=(800, 900),
                    title=["Epidemic Duration" "Maximum Number of Infected" "Susceptible Fraction Remaining"])
    savefig(combined_plot, "figures/combined_metrics_$(model.network_type)_mdeg_$(model.mean_degree).pdf")
    
    # Return comnined plot
    return combined_plot
end

"""
    compare_network_metrics(; network_types::Vector{Symbol}, mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1)

Compare key network metrics across different network types, visualized in a combined plot.

# Arguments
- `network_types::Vector{Symbol}`: Vector of network types to compare (e.g., [:random, :smallworld, :preferentialattachment])
- `mean_degree::Int`: The mean degree to use for all networks. Default is 4.
- `n_nodes::Int`: The number of nodes in each network. Default is 1000.
- `dispersion::Float64`: The dispersion parameter for the network. Default is 0.1.

# Returns
- A plot comparing network metrics across different network types

# Example
```julia
metrics_plot = compare_network_metrics(
    network_types=[:random, :smallworld, :preferentialattachment],
    mean_degree=4
)
```
"""
function compare_network_metrics(; network_types::Vector{Symbol}, mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1)
    # Store metrics for each network type
    all_metrics = Dict()
    
    # For each network type, generate a network and calculate its metrics
    for network_type in network_types
        println("Analyzing metrics for $(network_type) network...")
        model = initialize(; network_type, mean_degree, n_nodes, dispersion)
        analysis = analyze_graph(model.graph)
        
        # Extract metrics from the analysis results - handle the metric/value format
        summary_df = analysis["summary"]
        centrality_df = analysis["centrality"]
        
        # The summary dataframe has "metric" and "value" columns, so we extract values by filtering rows
        get_metric = function(metric_name)
            row = summary_df[summary_df.metric .== metric_name, :]
            return isempty(row) ? 0.0 : row.value[1]
        end
        
        # Store the summary metrics and centrality statistics
        all_metrics[network_type] = Dict(
            "density" => get_metric("network_density"),
            "clustering_coefficient" => get_metric("global_clustering_coefficient"),
            "assortativity" => get_metric("degree_assortativity"),
            "degree_centrality_mean" => mean(centrality_df.degree_centrality),
            "betweenness_centrality_mean" => mean(centrality_df.betweenness_centrality),
            "closeness_centrality_mean" => mean(centrality_df.closeness_centrality),
            "eigenvector_centrality_mean" => mean(centrality_df.eigenvector_centrality)
        )
        
        # Save the full results for this network type
        CSV.write("data/network_metrics_$(network_type)_mdeg_$(mean_degree)_nn_$(n_nodes).csv", summary_df)
        CSV.write("data/centrality_metrics_$(network_type)_mdeg_$(mean_degree)_nn_$(n_nodes).csv", centrality_df)
    end
    
    # Prepare data for plotting
    network_labels = [titlecase(String(nt)) for nt in network_types]
    
    # Create separate plots for different metric groups
    
    # 1. Structural metrics
    structure_metrics_plot = plot(
        ylabel="Value",
        title="Network Structure Metrics\n(mean degree: $(mean_degree), nodes: $(n_nodes))",
        legend=:topleft,
        size=(600, 400),
        margin=5mm,
        guidefontsize=9,
        titlefontsize=10,
        xtickfontsize=8
    )
    
    # Add density
    density_values = [all_metrics[nt]["density"] for nt in network_types]
    plot!(structure_metrics_plot, 1:length(network_types), density_values, label="Density", marker=:circle)
    
    # Add clustering coefficient
    clustering_values = [all_metrics[nt]["clustering_coefficient"] for nt in network_types]
    plot!(structure_metrics_plot, 1:length(network_types), clustering_values, label="Clustering Coefficient", marker=:square)
    
    # Add assortativity
    assortativity_values = [all_metrics[nt]["assortativity"] for nt in network_types]
    plot!(structure_metrics_plot, 1:length(network_types), assortativity_values, label="Assortativity", marker=:diamond)
    
    # Set x-axis ticks to network labels
    plot!(structure_metrics_plot, xticks=(1:length(network_types), network_labels))
    
    # 2. Centrality metrics
    centrality_metrics_plot = plot(
        ylabel="Mean Value",
        title="Mean Centrality Metrics\n(mean degree: $(mean_degree), nodes: $(n_nodes))",
        legend=:topleft,
        size=(600, 400),
        margin=5mm,
        guidefontsize=9,
        titlefontsize=10,
        xtickfontsize=8
    )
    
    # Add mean degree centrality
    degree_centrality_values = [all_metrics[nt]["degree_centrality_mean"] for nt in network_types]
    plot!(centrality_metrics_plot, 1:length(network_types), degree_centrality_values, label="Degree Centrality", marker=:circle)
    
    # Add mean betweenness centrality
    betweenness_centrality_values = [all_metrics[nt]["betweenness_centrality_mean"] for nt in network_types]
    plot!(centrality_metrics_plot, 1:length(network_types), betweenness_centrality_values, label="Betweenness Centrality", marker=:square)
    
    # Add mean closeness centrality
    closeness_centrality_values = [all_metrics[nt]["closeness_centrality_mean"] for nt in network_types]
    plot!(centrality_metrics_plot, 1:length(network_types), closeness_centrality_values, label="Closeness Centrality", marker=:diamond)
    
    # Add mean eigenvector centrality
    eigenvector_centrality_values = [all_metrics[nt]["eigenvector_centrality_mean"] for nt in network_types]
    plot!(centrality_metrics_plot, 1:length(network_types), eigenvector_centrality_values, label="Eigenvector Centrality", marker=:star5)
    
    # Set x-axis ticks to network labels
    plot!(centrality_metrics_plot, xticks=(1:length(network_types), network_labels))
    
    # Create a combined plot
    combined_metrics_plot = plot(structure_metrics_plot, centrality_metrics_plot, layout=(2,1), size=(800, 800))
    
    # Save the combined plot
    savefig(combined_metrics_plot, "figures/network_metrics_comparison_mdeg_$(mean_degree).pdf")
    
    return combined_metrics_plot
end