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
- `network_type::Symbol`: The type of network to use for the simulation. Possible values are `:random`, `:smallworld`, `:preferential`, `:configuration`, or `:proportionatemixing`.
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
    
    # Create a combined plot for side-by-side visualization
    combined_plot = plot(plotdynamics, plotdegdist, layout=(1,2), size=(1000, 400),
                       title=["Epidemic Dynamics ($(String(model.network_type)))" "Degree Distribution ($(String(model.network_type)))"])
    savefig(combined_plot, "figures/combined_plot_$(base_filename).pdf")
    
    # Return all three plots
    return plotdynamics, plotdegdist, combined_plot
end

"""
    create_comparison_box_plot(formatted_data, title, ylabel; filename=nothing)

Create a standardized box plot for comparison data.

# Arguments
- `formatted_data`: A DataFrame with columns `value` and `network_type`
- `title`: Plot title
- `ylabel`: Y-axis label
- `filename`: Optional filename to save the plot

# Returns
- A boxplot comparing the data across network types
"""
function create_comparison_box_plot(formatted_data, title, ylabel; filename=nothing, colors=nothing)
    # If colors are provided, use them for the boxplot
    p = boxplot(
        formatted_data.network_type, formatted_data.value,
        title=title,
        xlabel="Network Type",
        ylabel=ylabel,
        legend=false,
        size=(600, 300),
        margin=5mm,
        guidefontsize=9,
        titlefontsize=10,
        xtickfontsize=8,
        linewidth=1,
        left_margin=8mm,
        outliers=false
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
                          n_steps::Int=100, boxplot_colors=nothing)

Run simulations for multiple network types and generate comparison plots.
"""
function run_and_plot_comparison(; network_types::Vector{Symbol}, mean_degree::Int=4, 
                               n_nodes::Int=1000, dispersion::Float64=0.1, 
                               patient_zero::Symbol=:random, high_risk::Symbol=:random, 
                               fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, 
                               n_steps::Int=100, boxplot_colors=nothing)
    # Store results for each network type
    all_results = Dict()
    
    # Prepare empty DataFrames for each metric
    duration_data = DataFrame(value=Float64[], network_type=String[])
    max_infected_data = DataFrame(value=Float64[], network_type=String[])
    susceptible_remaining_data = DataFrame(value=Float64[], network_type=String[])
    
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
        println("Saving simulation results to data/simulation_results...")
        CSV.write("data/simulation_results_$(base_filename).csv", multiple_runs)
        println("Saving final results to output/final_results...")
        CSV.write("output/final_results_$(base_filename).csv", final_results)
        
        # Store results for box plots
        all_results[network_type] = final_results
        
        # Add data for each metric to the appropriate DataFrame with network type as a label
        network_label = titlecase(String(network_type))
        for value in final_results.first_to_last_infected
            push!(duration_data, (value, network_label))
        end
        
        for value in final_results.max_infected
            push!(max_infected_data, (value, network_label))
        end
        
        for value in final_results.susceptible_fraction_remaining
            push!(susceptible_remaining_data, (value, network_label))
        end
    end
    
    # Generate comparison plots using the pre-formatted data and colors
    duration_comparison = create_comparison_box_plot(
        duration_data,
        "Epidemic Duration Comparison\n(mean degree: $(mean_degree))",
        "Duration (steps)",
        filename="figures/duration_comparison_mdeg_$(mean_degree).pdf",
        colors=boxplot_colors
    )
    
    max_infected_comparison = create_comparison_box_plot(
        max_infected_data,
        "Maximum Infected Comparison\n(mean degree: $(mean_degree))",
        "Maximum Number of Infected",
        filename="figures/max_infected_comparison_mdeg_$(mean_degree).pdf",
        colors=boxplot_colors
    )
    
    sfr_comparison = create_comparison_box_plot(
        susceptible_remaining_data,
        "Susceptible Fraction Remaining Comparison\n(mean degree: $(mean_degree))",
        "Susceptible Fraction Remaining",
        filename="figures/sfr_comparison_mdeg_$(mean_degree).pdf",
        colors=boxplot_colors
    )
    
    # Create combined comparison plot
    combined_comparison = plot(
        duration_comparison, max_infected_comparison, sfr_comparison,
        layout=(1,3),
        size=(1200, 500),
        margin=3mm,
        title=["Epidemic Duration" "Maximum Infected" "Susceptible Fraction Remaining"],
        titlefontsize=10,
        plot_title="Epidemic Comparison (mean degree: $(mean_degree))", 
        plot_titlefontsize=12,
        left_margin=8mm,
    )
    println("Saving combined comparison plot to figures/combined_comparison_mdeg_$(mean_degree).pdf")
    savefig(combined_comparison, "figures/combined_comparison_mdeg_$(mean_degree).pdf")
    
    return combined_comparison
end

"""
    plot_network_metrics_comparison(metrics_data::Dict; save_path::String="")

Create a horizontal bar chart comparing network metrics across different network types.

# Arguments
- `metrics_data`: Dictionary with network types as keys and metrics as values
- `save_path`: Optional file path to save the plot

# Returns
- A horizontal bar chart comparing network metrics
"""
function plot_network_metrics_comparison(metrics_data::Dict; save_path::String="")
    # Get network types and prepare labels
    network_types = collect(keys(metrics_data))
    network_labels = Dict(
        :random => "Random",
        :smallworld => "Small-World", 
        :preferential => "Preferential Attachment",
        :configuration => "Configuration",
        :proportionate => "Proportionate Mixing"
    )
    
    # Define metrics to extract and their display names
    metrics = [
        "density", "clustering_coefficient", "assortativity",
        "degree_centrality", "betweenness_centrality", 
        "closeness_centrality", "eigenvector_centrality"
    ]
    
    metric_names = [
        "Density", "Clustering Coefficient", "Assortativity",
        "Degree Centrality", "Betweenness Centrality", 
        "Closeness Centrality", "Eigenvector Centrality"
    ]
    
    # Extract metrics for each network type
    data = Dict()
    
    for network_type in network_types
        analysis = metrics_data[network_type]
        
        data[network_type] = Dict(
            "density" => analysis["summary"][analysis["summary"].metric .== "Density", :value][1],
            "clustering_coefficient" => analysis["summary"][analysis["summary"].metric .== "Clustering Coefficient", :value][1],
            "assortativity" => analysis["summary"][analysis["summary"].metric .== "Assortativity", :value][1],
            "degree_centrality" => mean(analysis["centrality"].degree_centrality),
            "betweenness_centrality" => mean(analysis["centrality"].betweenness_centrality),
            "closeness_centrality" => mean(analysis["centrality"].closeness_centrality),
            "eigenvector_centrality" => mean(analysis["centrality"].eigenvector_centrality)
        )
    end
    
    # Prepare DataFrame for plotting
    df = DataFrame(metric = metric_names, position = 1:length(metric_names))
    
    # Add columns for each network type
    for network_type in network_types
        name = get(network_labels, network_type, String(network_type))
        df[!, name] = [data[network_type][metric] for metric in metrics]
    end
    
    # Set colors for network types
    colors = Dict(
        "Random" => :blue,
        "Small-World" => :orange,
        "Preferential Attachment" => :green
    )
    
    # Get network names
    network_names = setdiff(names(df), ["position", "metric"])
    
    # Display the data that will be plotted
    println("Plotting DataFrame:")
    println(df)
    
    # Create plot, using Vector{String} for yticks labels to ensure all display properly
    ytick_labels = Vector{String}(df.metric)
    
    # Properly normalize the data for better visualization
    # Scale the values to make them visible on the same plot
    scaled_df = deepcopy(df)
    for metric in metrics
        max_val = maximum([data[nt][metric] for nt in network_types])
        if max_val > 0
            idx = findfirst(==(metric_names[findfirst(==(metric), metrics)]), df.metric)
            if idx !== nothing
                # Apply logarithmic scaling to make small values more visible
                if max_val < 0.01
                    for network_name in network_names
                        scaled_df[idx, network_name] *= 100
                    end
                    ytick_labels[idx] = "$(ytick_labels[idx]) (×100)"
                end
            end
        end
    end
    
    p = StatsPlots.groupedbar(
        df.position,
        Matrix(df[:, network_names]),
        group=network_names,
        orientation=:horizontal,
        yticks=(df.position, ytick_labels),
        xlabel="Value",
        ylabel="Metric",
        title="Comparison of Network Metrics",
        legend=:bottomright,
        palette=[colors[n] for n in network_names],
        size=(900, 700),
        margin=10mm,
        guidefontsize=9,
        tickfontsize=8,
        titlefontsize=12
    )
    
    # Save if path provided
    if !isempty(save_path)
        savefig(p, save_path)
    end
    
    return p
end

"""
    compare_network_metrics(; network_types::Vector{Symbol}, mean_degree::Int=4, n_nodes::Int=1000)

Generate and compare metrics across different network types.

# Arguments
- `network_types`: List of network types to compare
- `mean_degree`: Mean degree for network generation
- `n_nodes`: Number of nodes in each network

# Returns
- A horizontal bar chart comparing network metrics
"""
function compare_network_metrics(; network_types::Vector{Symbol}, mean_degree::Int=4, n_nodes::Int=1000)
    # Generate and analyze each network type
    metrics_data = Dict()
    
    for nt in network_types
        println("Analyzing $(nt) network...")
        model = initialize(; network_type=nt, mean_degree=mean_degree, n_nodes=n_nodes)
        metrics_data[nt] = analyze_graph(model.graph)
    end
    
    # Generate and save the comparison plot
    filename = "figures/network_metrics_comparison_mdeg_$(mean_degree).pdf"
    plot = plot_network_metrics_comparison(metrics_data, save_path=filename)
    
    return plot
end

"""
    plot_centrality_comparison(;network_types=[:random, :smallworld, :preferential], mean_degree=4, n_nodes=1000, link_axes=false)

Plot boxplots comparing centrality measures across different network types.

# Arguments
- `network_types`: Vector of symbols representing the network types to compare (any number supported)
- `mean_degree`: Mean degree for network generation
- `n_nodes`: Number of nodes in each network
- `link_axes`: Boolean indicating whether to link y-axes across plots for easier comparison (default: false)

# Returns
- A combined plot showing boxplots of centrality measures for each network type

# Example
```julia
# Default visualization with independent y-axes for three network types
centrality_comparison = plot_centrality_comparison()

# With linked y-axes for direct comparison with two network types
centrality_comparison = plot_centrality_comparison(
    network_types=[:random, :preferential],
    link_axes=true
)

# With five different network types
centrality_comparison = plot_centrality_comparison(
    network_types=[:random, :smallworld, :preferential, :configuration, :proportionate]
)
```
"""
function plot_centrality_comparison(;network_types=[:random, :smallworld, :preferential], mean_degree=4, n_nodes=1000, link_axes=false)
    # Initialize empty DataFrames to store the centrality data
    centrality_data = Dict()
    
    # Generate and analyze each network type
    for nt in network_types
        model = initialize(; network_type=nt, mean_degree=mean_degree, n_nodes=n_nodes)
        analysis = analyze_graph(model.graph)
        centrality_data[nt] = analysis["centrality"]
    end
    
    # Determine the optimal plot layout based on the number of network types
    n_types = length(network_types)
    
    if n_types == 1
        # For a single network type, use 1x1
        plot_layout = (1, 1)
        plot_width = 600
    elseif n_types == 2
        # For two network types, use 1x2
        plot_layout = (1, 2)
        plot_width = 1000
    elseif n_types <= 4
        # For 3-4 network types, use 1xN
        plot_layout = (1, n_types)
        plot_width = min(1600, 500 * n_types)
    else
        # For more than 4 network types, use a more compact grid layout
        n_cols = ceil(Int, sqrt(n_types))
        n_rows = ceil(Int, n_types / n_cols)
        plot_layout = (n_rows, n_cols)
        plot_width = min(1800, 450 * n_cols)
    end
    
    # Fixed plot height per row
    plot_height_per_row = 450
    plot_height = plot_height_per_row * first(plot_layout)
    
    # Create plots for each centrality measure
    plots = []
    
    # Calculate y-axis limits for each network type independently
    # This allows better visualization of the specific distributions
    y_limits = Dict()
    
    for nt in network_types
        df = centrality_data[nt]
        
        # Calculate quartiles for each centrality measure
        q1_degree = quantile(df.degree_centrality, 0.25)
        q3_degree = quantile(df.degree_centrality, 0.75)
        iqr_degree = q3_degree - q1_degree
        upper_whisker_degree = min(maximum(df.degree_centrality), q3_degree + 1.5 * iqr_degree)
        
        q1_betweenness = quantile(df.betweenness_centrality, 0.25)
        q3_betweenness = quantile(df.betweenness_centrality, 0.75)
        iqr_betweenness = q3_betweenness - q1_betweenness
        upper_whisker_betweenness = min(maximum(df.betweenness_centrality), q3_betweenness + 1.5 * iqr_betweenness)
        
        q1_closeness = quantile(df.closeness_centrality, 0.25)
        q3_closeness = quantile(df.closeness_centrality, 0.75)
        iqr_closeness = q3_closeness - q1_closeness
        upper_whisker_closeness = min(maximum(df.closeness_centrality), q3_closeness + 1.5 * iqr_closeness)
        
        q1_eigenvector = quantile(df.eigenvector_centrality, 0.25)
        q3_eigenvector = quantile(df.eigenvector_centrality, 0.75)
        iqr_eigenvector = q3_eigenvector - q1_eigenvector
        upper_whisker_eigenvector = min(maximum(df.eigenvector_centrality), q3_eigenvector + 1.5 * iqr_eigenvector)
        
        # Use upper whiskers (excluding outliers) with a small buffer for y-axis limits
        degree_max = upper_whisker_degree * 1.15
        betweenness_max = upper_whisker_betweenness * 1.15
        closeness_max = upper_whisker_closeness * 1.15
        eigenvector_max = upper_whisker_eigenvector * 1.15
        
        # Store the maximum overall for this network type
        y_limits[nt] = max(degree_max, betweenness_max, closeness_max, eigenvector_max)
    end
    
    # If link_axes is true, find the global maximum across all networks
    global_y_max = link_axes ? maximum(values(y_limits)) : 0
    
    # Create a DataFrame for each network type to use with StatsPlots groupedboxplot
    for (i, nt) in enumerate(network_types)
        df = centrality_data[nt]
        network_name = titlecase(String(nt))
        
        # Create a DataFrame for boxplot data
        boxplot_data = DataFrame(
            centrality_type = vcat(
                fill("Degree", nrow(df)),
                fill("Betweenness", nrow(df)),
                fill("Closeness", nrow(df)),
                fill("Eigenvector", nrow(df))
            ),
            value = vcat(
                df.degree_centrality,
                df.betweenness_centrality,
                df.closeness_centrality,
                df.eigenvector_centrality
            )
        )
        
        # Convert centrality_type to a categorical array with ordered levels
        boxplot_data.centrality_type = CategoricalArray(
            boxplot_data.centrality_type, 
            ordered=true, 
            levels=["Degree", "Betweenness", "Closeness", "Eigenvector"]
        )
        
        # Set y-axis limits based on link_axes parameter
        y_max = link_axes ? global_y_max : y_limits[nt]
        
        # Create a more precisely aligned boxplot using @df macro
        p = @df boxplot_data boxplot(
            :centrality_type, 
            :value,
            title=network_name,
            legend=false,
            outliers=false,  # Hide outliers to improve readability
            marker=(0.5, :circle, 0.3),
            alpha=0.7,       
            guidefontsize=9,
            titlefontsize=10,
            widen=true,      
            bar_width=0.7,   
            xticks=([1,2,3,4], ["Degree", "Betweenness", "Closeness", "Eigenvector"]),
            xrotation=30,    
            ylims=(0, y_max),
            ylabel="Centrality Value",
            dpi=300,
            margin=8mm,      # Add individual plot margins for better spacing
            bottom_margin=10mm,  # Add extra space at the bottom for x-axis labels
            left_margin=8mm      # Add extra space for y-axis labels
        )
        
        push!(plots, p)
    end
    
    # Handle case where we might have more layout slots than plots
    if length(plots) < prod(plot_layout)
        # Fill remaining slots with empty plots
        for i in length(plots)+1:prod(plot_layout)
            push!(plots, plot(framestyle=:none, grid=false, showaxis=false))
        end
    end
    
    # Combine plots in a single figure with improved margins
    combined_plot = plot(plots..., 
                      layout=plot_layout, 
                      size=(plot_width, plot_height),
                      bottom_margin=12mm, # Increase bottom margin for x-axis labels
                      top_margin=10mm,    # Add more top margin for titles
                      xtickfontsize=9,
                      link=link_axes ? :y : :none,
                      title_position=:center)
    
    # Save the figure
    net_types_str = join([String(nt) for nt in network_types], "_")
    println("Saving combined centrality comparison plot to figures/centrality_comparison_$(net_types_str)_mdeg_$(mean_degree).pdf")
    savefig(combined_plot, "figures/centrality_comparison_$(net_types_str)_mdeg_$(mean_degree).pdf")
    
    return combined_plot
end