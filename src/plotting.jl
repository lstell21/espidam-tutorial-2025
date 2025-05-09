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
    # Extract degrees and counts from different possible input formats
    if isa(degree_distribution, Dict)
        degrees = collect(keys(degree_distribution))
        counts = collect(values(degree_distribution))
    elseif isa(degree_distribution, DataFrame) && haskey(degree_distribution, :degree_distribution)
        dd = degree_distribution.degree_distribution[1]
        degrees = collect(keys(dd))
        counts = collect(values(dd))
    elseif isa(degree_distribution, DataFrame) && haskey(degree_distribution, :degree) && haskey(degree_distribution, :count)
        degrees = degree_distribution.degree
        counts = degree_distribution.count
    else
        error("Unsupported degree_distribution format")
    end
    
    # Create and return the plot
    return plot(degrees, counts, seriestype=:bar, 
                title="Degree Distribution ($(String(network_type)))", 
                xlabel="Degree", ylabel="Count", legend=:none)
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
    # Extract the data and create time vector
    susceptible = mdf[!, :susceptible_count]
    infected = mdf[!, :infected_count]
    recovered = mdf[!, :recovered_count]
    time = 1:length(susceptible)

    # Create and return the plot
    p = plot(time, susceptible, label="Susceptible", legend=:topright, 
             xlabel="Time", ylabel="Count", linewidth=2,
             title="Epidemic Dynamics ($(String(network_type)))")
    plot!(p, time, infected, label="Infected", linewidth=2)
    plot!(p, time, recovered, label="Recovered", linewidth=2)
    
    return p
end

"""
    create_base_filename(model)

Create a consistent base filename from model parameters for saving files.
"""
function create_base_filename(model)
    return "$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob)"
end

"""
    plot_single_run(; network_type::Symbol, mean_degree::Int=4, n_nodes::Int=1000, 
                   dispersion::Float64=0.1, patient_zero::Symbol=:random, 
                   high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, 
                   trans_prob::Float64=0.1, n_steps::Int=100)

Plot a single run of an epidemic simulation.

# Arguments
- `network_type::Symbol`: The type of network to use for the simulation. Possible values are `:random`, `:smallworld`, `:preferentialattachment`, `:configuration`, or `:proportionatemixing`.
- `mean_degree::Int`: The mean degree of the network. Default is 4.
- `n_nodes::Int`: The number of nodes in the network. Default is 1000.
- `dispersion::Float64`: The dispersion parameter for the network. Default is 0.1.
- `patient_zero::Symbol`: The type of patient zero to use for the simulation. Default is `:random`.
- `high_risk::Symbol`: How high-risk individuals are distributed. Default is `:random`.
- `fraction_high_risk::Float64`: The fraction of high-risk individuals in the population. Default is 0.1.
- `trans_prob::Float64`: The transmission probability. Default is 0.1.
- `n_steps::Int`: The number of simulation steps to run. Default is 100.

# Returns
- `plotdynamics`: A plot of the epidemic trajectories.
- `plotdegdist`: A plot of the degree distribution.
- `combined_plot`: A combined plot with both plots side-by-side.

# Example
```julia
dynamics_plot, degdist_plot, combined_plot = plot_single_run(network_type=:random, mean_degree=2)
```
"""
function plot_single_run(; network_type::Symbol, mean_degree::Int=4, n_nodes::Int=1000, 
                       dispersion::Float64=0.1, patient_zero::Symbol=:random, 
                       high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, 
                       trans_prob::Float64=0.1, n_steps::Int=100)
    # Initialize model and run simulation
    model = initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero, 
                     high_risk, fraction_high_risk, trans_prob)
    
    # Define adata and mdata locally to avoid relying on global variables
    adata = [:status]
    mdata = [:susceptible_count, :infected_count, :recovered_count]
    
    # Run the simulation
    _, mdf = run!(model, n_steps; adata, mdata)
    
    # Get graph analysis results
    graph_analysis = analyze_graph(model.graph)
    
    # Create base filename for consistency
    base_filename = create_base_filename(model)
    
    # Create plots
    plotdynamics = plot_epidemic_trajectories(mdf, model.network_type)
    plotdegdist = plot_degree_distribution(graph_analysis["degree_distribution"]; network_type=model.network_type)
    
    # Save plots
    savefig(plotdynamics, "figures/plotdynamics_$(base_filename).pdf")
    savefig(plotdegdist, "figures/plotdegdist_$(base_filename).pdf")
    
    # Create a combined plot for side-by-side visualization
    combined_plot = plot(plotdynamics, plotdegdist, layout=(1,2), size=(1000, 400),
                       title=["Epidemic Dynamics ($(String(model.network_type)))" "Degree Distribution ($(String(model.network_type)))"])
    savefig(combined_plot, "figures/combined_plot_$(base_filename).pdf")
    
    # Return all three plots
    return plotdynamics, plotdegdist, combined_plot
end

"""
    create_comparison_bar_plot(means, stds, labels, title, ylabel; filename=nothing)

Create a standardized bar plot for comparison data.
"""
function create_comparison_bar_plot(means, stds, labels, title, ylabel; filename=nothing)
    p = bar(
        labels, means, yerror=stds,
        title=title,
        xlabel="Network Type",
        ylabel=ylabel,
        legend=false,
        size=(600, 300),
        margin=5mm,
        guidefontsize=9,
        titlefontsize=10,
        xtickfontsize=8,
        bar_width=0.4,
        linewidth=1,
        left_margin=8mm
    )
    
    if !isnothing(filename)
        savefig(p, filename)
    end
    
    return p
end

"""
    run_and_plot_comparison(; network_types::Vector{Symbol}, mean_degree::Int=4, 
                          n_nodes::Int=1000, dispersion::Float64=0.1, 
                          patient_zero::Symbol=:random, high_risk::Symbol=:random, 
                          fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, 
                          n_steps::Int=100)

Run simulations for multiple network types and generate comparison plots.
"""
function run_and_plot_comparison(; network_types::Vector{Symbol}, mean_degree::Int=4, 
                               n_nodes::Int=1000, dispersion::Float64=0.1, 
                               patient_zero::Symbol=:random, high_risk::Symbol=:random, 
                               fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, 
                               n_steps::Int=100)
    # Store results for each network type
    all_results = Dict()
    
    # Run simulations for each network type
    for network_type in network_types
        println("Running simulations for $(network_type) network...")
        model = initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero, 
                          high_risk, fraction_high_risk, trans_prob)
        
        # Run simulations
        multiple_runs = run_simulations(; network_type, mean_degree, n_nodes, dispersion, 
                                       patient_zero, high_risk, fraction_high_risk, 
                                       trans_prob, n_steps)
        
        # Process results
        grouped_data = groupby(multiple_runs, [:seed])
        final_results = combine(grouped_data, :infected_count => argmin => :first_to_last_infected, 
                               :infected_count => maximum => :max_infected)
        last_rows = combine(grouped_data, names(multiple_runs) .=> last)
        final_results[!, :susceptible_fraction_remaining] = last_rows.susceptible_count_last ./ 
            (last_rows.susceptible_count_last + last_rows.infected_count_last + last_rows.recovered_count_last)
        
        # Save results
        base_filename = "$(network_type)_mdeg_$(mean_degree)_nn_$(n_nodes)_disp_$(dispersion)_pat0_$(patient_zero)_hirisk_$(high_risk)_hr_frac_$(fraction_high_risk)_trans_$(trans_prob)"
        CSV.write("data/simulation_results_$(base_filename).csv", multiple_runs)
        CSV.write("output/final_results_$(base_filename).csv", final_results)
        
        # Calculate summary statistics
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
    
    # Extract metrics for plotting
    duration_means = [all_results[nt].duration_mean for nt in network_types]
    duration_stds = [all_results[nt].duration_std for nt in network_types]
    max_infected_means = [all_results[nt].max_infected_mean for nt in network_types]
    max_infected_stds = [all_results[nt].max_infected_std for nt in network_types]
    sfr_means = [all_results[nt].sfr_mean for nt in network_types]
    sfr_stds = [all_results[nt].sfr_std for nt in network_types]
    
    # Generate comparison plots
    duration_comparison = create_comparison_bar_plot(
        duration_means, duration_stds, network_labels, 
        "Epidemic Duration Comparison\n(mean degree: $(mean_degree))",
        "Duration (steps)",
        filename="figures/duration_comparison_mdeg_$(mean_degree).pdf"
    )
    
    max_infected_comparison = create_comparison_bar_plot(
        max_infected_means, max_infected_stds, network_labels, 
        "Maximum Infected Comparison\n(mean degree: $(mean_degree))",
        "Maximum Number of Infected",
        filename="figures/max_infected_comparison_mdeg_$(mean_degree).pdf"
    )
    
    sfr_comparison = create_comparison_bar_plot(
        sfr_means, sfr_stds, network_labels, 
        "Susceptible Fraction Remaining Comparison\n(mean degree: $(mean_degree))",
        "Susceptible Fraction Remaining",
        filename="figures/sfr_comparison_mdeg_$(mean_degree).pdf"
    )
    
    # Create combined comparison plot
    combined_comparison = plot(
        duration_comparison, max_infected_comparison, sfr_comparison,
        layout=(3,1),
        size=(700, 700),
        margin=3mm,
        title=["Epidemic Duration" "Maximum Infected" "Susceptible Fraction Remaining"],
        titlefontsize=10,
        plot_title="Network Comparison (mean degree: $(mean_degree))", 
        plot_titlefontsize=12,
        left_margin=8mm
    )
    savefig(combined_comparison, "figures/combined_comparison_mdeg_$(mean_degree).pdf")
    
    return combined_comparison
end

"""
    run_and_plot(; network_type::Symbol, mean_degree::Int=4, n_nodes::Int=1000, 
                dispersion::Float64=0.1, patient_zero::Symbol=:random, 
                high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, 
                trans_prob::Float64=0.1, n_steps::Int=100)

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
    final_results = combine(grouped_data, :infected_count => argmin => :first_to_last_infected, 
                           :infected_count => maximum => :max_infected)
    last_rows = combine(grouped_data, names(multiple_runs) .=> last)
    final_results[!, :susceptible_fraction_remaining] = last_rows.susceptible_count_last ./ 
        (last_rows.susceptible_count_last + last_rows.infected_count_last + last_rows.recovered_count_last)
    
    # Save results
    base_filename = create_base_filename(model)
    CSV.write("data/simulation_results_$(base_filename).csv", multiple_runs)
    CSV.write("output/final_results_$(base_filename).csv", final_results)
    
    # Create descriptive title
    network_desc = "$(titlecase(String(network_type))) Network (mean degree: $(mean_degree))"
    
    # Create plots using a common function to reduce redundancy
    plot_bar_data = function(data, ylabel, subtitle)
        plot(data, 
             seriestype=:bar, 
             xlabel="Simulation Run", 
             ylabel=ylabel,
             legend=:none,
             title="$(subtitle)\n$(network_desc)",
             guidefontsize=10,
             titlefontsize=12,
             size=(700, 500),
             margin=10mm)
    end
    
    # Generate the three plots
    duration_plot = plot_bar_data(final_results.first_to_last_infected, 
                              "Duration of Epidemic (steps)", "Epidemic Duration")
    savefig(duration_plot, "figures/duration_plot_$(model.network_type)_mdeg_$(model.mean_degree).pdf")
    
    max_infected_plot = plot_bar_data(final_results.max_infected, 
                                  "Maximum Number of Infected", "Maximum Number of Infected")
    savefig(max_infected_plot, "figures/max_infected_plot_$(model.network_type)_mdeg_$(model.mean_degree).pdf")
    
    sfr_plot = plot_bar_data(final_results.susceptible_fraction_remaining, 
                         "Susceptible Fraction Remaining", "Susceptible Fraction Remaining")
    savefig(sfr_plot, "figures/sfr_plot_$(model.network_type)_mdeg_$(model.mean_degree).pdf")
    
    # Create combined plot
    combined_plot = plot(duration_plot, max_infected_plot, sfr_plot, 
                    layout=(3,1), 
                    size=(800, 900),
                    title=["Epidemic Duration" "Maximum Number of Infected" "Susceptible Fraction Remaining"])
    savefig(combined_plot, "figures/combined_metrics_$(model.network_type)_mdeg_$(model.mean_degree).pdf")
    
    return combined_plot
end

"""
    compare_network_metrics(; network_types::Vector{Symbol}, mean_degree::Int=4, n_nodes::Int=1000, dispersion::Float64=0.1, use_log_scale::Bool=false, separate_assortativity::Bool=false, special_handle_clustering::Bool=false)

Compare key network metrics across different network types, visualized in a combined plot.

# Arguments
- `network_types::Vector{Symbol}`: Vector of network types to compare (e.g., [:random, :smallworld, :preferentialattachment])
- `mean_degree::Int`: The mean degree to use for all networks. Default is 4.
- `n_nodes::Int`: The number of nodes in each network. Default is 1000.
- `dispersion::Float64`: The dispersion parameter for the network. Default is 0.1.
- `use_log_scale::Bool`: Whether to use a log scale for the y-axis. Default is false.
- `separate_assortativity::Bool`: Whether to separate assortativity into its own panel. Default is false.
- `special_handle_clustering::Bool`: Whether to apply special handling for clustering coefficient visualization. Default is false.

# Returns
- A plot comparing network metrics across different network types

# Example
```julia
metrics_plot = compare_network_metrics(
    network_types=[:random, :smallworld, :preferentialattachment],
    mean_degree=4,
    use_log_scale=true,
    separate_assortativity=true,
    special_handle_clustering=true
)
```
"""
function compare_network_metrics(; network_types::Vector{Symbol}, mean_degree::Int=4, 
                               n_nodes::Int=1000, dispersion::Float64=0.1,
                               use_log_scale::Bool=false, separate_assortativity::Bool=false)
    # Store metrics for each network type
    all_metrics = Dict()
    
    # Helper function to extract metric from summary dataframe
    get_metric = function(summary_df, metric_name)
        row = summary_df[summary_df.metric .== metric_name, :]
        return isempty(row) ? 0.0 : row.value[1]
    end
    
    # Analyze each network type
    for network_type in network_types
        println("Analyzing metrics for $(network_type) network...")
        model = initialize(; network_type, mean_degree, n_nodes, dispersion)
        analysis = analyze_graph(model.graph)
        
        # Extract metrics
        summary_df = analysis["summary"]
        centrality_df = analysis["centrality"]
        
        # Store metrics
        all_metrics[network_type] = Dict(
            "density" => get_metric(summary_df, "Density"),
            "clustering_coefficient" => get_metric(summary_df, "Clustering Coefficient"),
            "assortativity" => get_metric(summary_df, "Assortativity"),
            "degree_centrality_mean" => mean(centrality_df.degree_centrality),
            "betweenness_centrality_mean" => mean(centrality_df.betweenness_centrality),
            "closeness_centrality_mean" => mean(centrality_df.closeness_centrality),
            "eigenvector_centrality_mean" => mean(centrality_df.eigenvector_centrality)
        )
        
        # Save results
        CSV.write("data/network_metrics_$(network_type)_mdeg_$(mean_degree)_nn_$(n_nodes).csv", summary_df)
        CSV.write("data/centrality_metrics_$(network_type)_mdeg_$(mean_degree)_nn_$(n_nodes).csv", centrality_df)
    end
    
    # Prepare data for plotting
    network_labels = [titlecase(String(nt)) for nt in network_types]
    
    # Extract all metrics
    density_values = [all_metrics[nt]["density"] for nt in network_types]
    clustering_values = [all_metrics[nt]["clustering_coefficient"] for nt in network_types]
    assortativity_values = [all_metrics[nt]["assortativity"] for nt in network_types]
    degree_centrality_values = [all_metrics[nt]["degree_centrality_mean"] for nt in network_types]
    betweenness_centrality_values = [all_metrics[nt]["betweenness_centrality_mean"] for nt in network_types]
    closeness_centrality_values = [all_metrics[nt]["closeness_centrality_mean"] for nt in network_types]
    eigenvector_centrality_values = [all_metrics[nt]["eigenvector_centrality_mean"] for nt in network_types]
    
    if separate_assortativity
        # Create a multi-panel plot with separate handling for assortativity
        
        # Panel 1: Positive-valued metrics with log scale (if requested)
        clustering_label = "Clustering Coefficient"
        
        positive_metrics_plot = plot(
            title="Network Metrics (Positive-valued)",
            ylabel=use_log_scale ? "Value (Log Scale)" : "Value",
            legend=:topleft,
            yscale=use_log_scale ? :log10 : :identity,
            size=(600, 400),
            margin=5mm
        )
        
        # Make sure all values for log scale are positive
        if use_log_scale
            density_values = [max(v, 1e-10) for v in density_values]
            clustering_values = [max(v, 1e-10) for v in clustering_values]
            degree_centrality_values = [max(v, 1e-10) for v in degree_centrality_values]
            betweenness_centrality_values = [max(v, 1e-10) for v in betweenness_centrality_values]
            closeness_centrality_values = [max(v, 1e-10) for v in closeness_centrality_values]
            eigenvector_centrality_values = [max(v, 1e-10) for v in eigenvector_centrality_values]
        end
        
        plot!(positive_metrics_plot, 1:length(network_types), density_values, 
              label="Density", marker=:circle, markersize=6, linewidth=2)
        plot!(positive_metrics_plot, 1:length(network_types), clustering_values, 
              label=clustering_label, marker=:square, markersize=6, linewidth=2)
        plot!(positive_metrics_plot, 1:length(network_types), degree_centrality_values, 
              label="Degree Centrality", marker=:star5, markersize=6, linewidth=2)
        plot!(positive_metrics_plot, 1:length(network_types), betweenness_centrality_values, 
              label="Betweenness Centrality", marker=:utriangle, markersize=6, linewidth=2)
        plot!(positive_metrics_plot, 1:length(network_types), closeness_centrality_values, 
              label="Closeness Centrality", marker=:dtriangle, markersize=6, linewidth=2)
        plot!(positive_metrics_plot, 1:length(network_types), eigenvector_centrality_values, 
              label="Eigenvector Centrality", marker=:hexagon, markersize=6, linewidth=2)
        plot!(positive_metrics_plot, xticks=(1:length(network_types), network_labels))
        
        # Panel 2: Assortativity (which can be negative)
        assortativity_plot = plot(
            title="Degree Assortativity",
            ylabel="Assortativity Value",
            legend=false,
            size=(600, 300),
            margin=5mm
        )
        
        # Use bar plot for assortativity to clearly show positive/negative values
        bar!(assortativity_plot, 1:length(network_types), assortativity_values, 
             fillcolor=:lightblue, linecolor=:blue, alpha=0.7)
        plot!(assortativity_plot, xticks=(1:length(network_types), network_labels))
        hline!(assortativity_plot, [0], linestyle=:dash, color=:gray, label="")
        
        # Combine the plots
        combined_metrics_plot = plot(
            positive_metrics_plot, assortativity_plot,
            layout=(2, 1),
            size=(800, 700),
            plot_title="Network Metrics Comparison (mean degree: $(mean_degree))",
            plot_titlefontsize=14
        )
        
        # Save plot
        log_suffix = use_log_scale ? "_log" : ""
        savefig(combined_metrics_plot, "figures/network_metrics_comparison$(log_suffix)_mdeg_$(mean_degree).pdf")
        
        return combined_metrics_plot
    else
        # Original single-plot approach with all metrics together
        # Adjustments for clustering if needed
        clustering_label = "Clustering Coefficient"
        
        # Create a unified plot for all metrics
        scale_type = use_log_scale ? :log10 : :identity
        scale_label = use_log_scale ? "(Log Scale)" : ""
        
        combined_metrics_plot = plot(
            ylabel="Value $scale_label",
            title="Network Metrics Comparison $scale_label\n(mean degree: $(mean_degree), nodes: $(n_nodes))",
            legend=:topleft,
            size=(800, 500),
            margin=8mm,
            guidefontsize=9,
            titlefontsize=12,
            xtickfontsize=10,
            linewidth=2,
            yscale=scale_type
        )
        
        # Ensure positive values for log scale if needed
        if use_log_scale
            # Handle each metric separately
            density_values = [max(v, 1e-10) for v in density_values]
            clustering_values = [max(v, 1e-10) for v in clustering_values]
            # Skip assortativity for log scale as it can be negative
            degree_centrality_values = [max(v, 1e-10) for v in degree_centrality_values]
            betweenness_centrality_values = [max(v, 1e-10) for v in betweenness_centrality_values]
            closeness_centrality_values = [max(v, 1e-10) for v in closeness_centrality_values]
            eigenvector_centrality_values = [max(v, 1e-10) for v in eigenvector_centrality_values]
        end
        
        # Add metrics to plot (excluding assortativity for log scale)
        plot!(combined_metrics_plot, 1:length(network_types), density_values, 
            label="Density", marker=:circle, markersize=6)
        plot!(combined_metrics_plot, 1:length(network_types), clustering_values, 
            label=clustering_label, marker=:square, markersize=6)
        
        if !use_log_scale || !any(assortativity_values .< 0)
            plot!(combined_metrics_plot, 1:length(network_types), assortativity_values, 
                label="Assortativity", marker=:diamond, markersize=6)
        end
        
        plot!(combined_metrics_plot, 1:length(network_types), degree_centrality_values, 
            label="Degree Centrality", marker=:star5, markersize=6)
        plot!(combined_metrics_plot, 1:length(network_types), betweenness_centrality_values, 
            label="Betweenness Centrality", marker=:utriangle, markersize=6)
        plot!(combined_metrics_plot, 1:length(network_types), closeness_centrality_values, 
            label="Closeness Centrality", marker=:dtriangle, markersize=6)
        plot!(combined_metrics_plot, 1:length(network_types), eigenvector_centrality_values, 
            label="Eigenvector Centrality", marker=:hexagon, markersize=6)
        
        # Set x-axis labels to network types
        plot!(combined_metrics_plot, xticks=(1:length(network_types), network_labels))
        
        # Save plot
        log_suffix = use_log_scale ? "_log" : ""
        savefig(combined_metrics_plot, "figures/network_metrics_comparison$(log_suffix)_mdeg_$(mean_degree).pdf")
        
        return combined_metrics_plot
    end
end