using CategoricalArrays

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
        CSV.write("data/simulation_results_$(base_filename).csv", multiple_runs)
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
        plot_title="Network Comparison (mean degree: $(mean_degree))", 
        plot_titlefontsize=12,
        left_margin=8mm,
    )
    savefig(combined_comparison, "figures/combined_comparison_mdeg_$(mean_degree).pdf")
    
    return combined_comparison
end

"""
    plot_network_metrics_comparison(metrics_data::Dict; save_path::Union{String, Nothing}=nothing)

Create a horizontal bar chart comparing various network metrics across different network types.

# Arguments
- `metrics_data::Dict`: Dictionary with network types as keys and their metrics as values
- `save_path::Union{String, Nothing}=nothing`: Optional file path to save the plot

# Returns
- A horizontal bar chart comparing network metrics across different network types

# Example
```julia
# Create networks and analyze them
network_types = [:random, :smallworld, :preferential]
metrics_data = Dict()
for nt in network_types
    model = initialize(; network_type=nt, mean_degree=4)
    metrics_data[nt] = analyze_graph(model.graph)
end
# Plot comparison
plot_network_metrics_comparison(metrics_data)
```
"""
function plot_network_metrics_comparison(metrics_data::Dict; save_path::Union{String, Nothing}=nothing)
    # Network types to include in the comparison
    network_types = collect(keys(metrics_data))
    network_labels = Dict(
        :random => "random",
        :smallworld => "smallworld", 
        :preferentialattachment => "preferentialattachment"
    )
    
    # Define metrics to extract and their display names
    metrics = [
        "density",
        "clustering_coefficient",
        "assortativity",
        "degree_centrality",
        "betweenness_centrality",
        "closeness_centrality",
        "eigenvector_centrality"
    ]
    
    metric_display_names = [
        "Density",
        "Clustering Coefficient",
        "Assortativity",
        "Degree Centrality",
        "Betweenness Centrality",
        "Closeness Centrality", 
        "Eigenvector Centrality"
    ]
    
    # Create empty dictionary to store extracted metrics
    extracted_metrics = Dict()
    
    # For each network type, extract the relevant metrics
    for network_type in network_types
        analysis = metrics_data[network_type]
        
        # Get values from analysis
        centrality_df = analysis["centrality"]
        summary_df = analysis["summary"]
        
        # Extract values
        extracted_metrics[network_type] = Dict(
            "density" => summary_df[summary_df.metric .== "Density", :value][1],
            "clustering_coefficient" => summary_df[summary_df.metric .== "Clustering Coefficient", :value][1],
            "assortativity" => summary_df[summary_df.metric .== "Assortativity", :value][1],
            "degree_centrality" => mean(centrality_df.degree_centrality),
            "betweenness_centrality" => mean(centrality_df.betweenness_centrality),
            "closeness_centrality" => mean(centrality_df.closeness_centrality),
            "eigenvector_centrality" => mean(centrality_df.eigenvector_centrality)
        )
    end
    
    # Prepare data for plotting - create a DataFrame for grouped bar plotting
    plot_data = DataFrame()
    
    # Create a position column for the bars
    plot_data.position = 1:length(metric_display_names)
    plot_data.metric = reverse(metric_display_names)
    
    # Add a column for each network type
    for network_type in network_types
        network_name = get(network_labels, network_type, String(network_type))
        values = Float64[]
        
        # Get values in the same order as metric_display_names (reversed)
        for metric in reverse(metrics)
            push!(values, extracted_metrics[network_type][metric])
        end
        
        # Add the column to the DataFrame
        plot_data[!, network_name] = values
    end
    
    # Define colors for each network type
    colors = Dict(
        "random" => :blue,
        "smallworld" => :orange,
        "preferentialattachment" => :green
    )
    
    # Get network names from the DataFrame columns (excluding position and metric)
    network_names = setdiff(names(plot_data), ["position", "metric"])
    
    # Create horizontal bar chart - Using direct array access instead of col() function
    p = groupedbar(
        plot_data.position,
        Matrix(plot_data[:, network_names]),  # Directly accessing columns as a matrix
        group=network_names,
        orientation=:horizontal,
        yticks=(plot_data.position, plot_data.metric),  # y-axis ticks with metric names
        xlabel="Value",
        ylabel="Metric",
        title="Comparison of Graph Measures",
        legend=:bottomright,
        palette=[colors[n] for n in network_names],
        size=(800, 600),
        margin=10mm,
        left_margin=18mm,
        grid=true,
        framestyle=:box,
        linewidth=1
    )
    
    # Save if path provided
    if !isnothing(save_path)
        savefig(p, save_path)
    end
    
    return p
end

"""
    compare_network_metrics(; network_types::Vector{Symbol}=[:random, :smallworld, :preferentialattachment], 
                          mean_degree::Int=4, n_nodes::Int=1000)

Generate, analyze, and compare metrics across different network types.

# Arguments
- `network_types::Vector{Symbol}`: List of network types to compare
- `mean_degree::Int`: Mean degree for network generation
- `n_nodes::Int`: Number of nodes in each network

# Returns
- A horizontal bar chart comparing network metrics across different network types

# Example
```julia
compare_network_metrics(mean_degree=6)
```
"""
function compare_network_metrics(; network_types::Vector{Symbol}=[:random, :smallworld, :preferentialattachment], 
                               mean_degree::Int=4, n_nodes::Int=1000)
    
    metrics_data = Dict()
    
    # Generate and analyze each network type
    for nt in network_types
        println("Analyzing $(nt) network...")
        model = initialize(; network_type=nt, mean_degree=mean_degree, n_nodes=n_nodes)
        metrics_data[nt] = analyze_graph(model.graph)
    end
    
    # Generate the comparison plot
    plot = plot_network_metrics_comparison(metrics_data, 
                                         save_path="figures/network_metrics_comparison_mdeg_$(mean_degree).pdf")
    
    return plot
end

"""
    plot_centrality_comparison(;network_types=[:random, :smallworld, :preferentialattachment], mean_degree=4, n_nodes=1000)

Plot boxplots comparing centrality measures across different network types.

# Arguments
- `network_types`: Vector of symbols representing the network types to compare
- `mean_degree`: Mean degree for network generation
- `n_nodes`: Number of nodes in each network

# Returns
- A combined plot showing boxplots of centrality measures for each network type

# Example
```julia
centrality_comparison = plot_centrality_comparison()
display(centrality_comparison)
```
"""
function plot_centrality_comparison(;network_types=[:random, :smallworld, :preferentialattachment], mean_degree=4, n_nodes=1000)
    # Initialize empty DataFrames to store the centrality data
    centrality_data = Dict()
    
    # Generate and analyze each network type
    for nt in network_types
        model = initialize(; network_type=nt, mean_degree=mean_degree, n_nodes=n_nodes)
        analysis = analyze_graph(model.graph)
        centrality_data[nt] = analysis["centrality"]
    end
    
    # Setup figure layout - use a simple grid layout instead of percentage-based layout
    plot_layout = (1, 3)
    
    # Colors for the different network types
    colors = Dict(:random => :blue, :smallworld => :gray, :preferentialattachment => :orange)
    
    # Create plots for each centrality measure
    plots = []
    
    # Calculate y-axis limits for each network type independently
    # This allows better visualization of the specific distributions
    y_limits = Dict()
    
    for nt in network_types
        df = centrality_data[nt]
        
        # Find the maximum values for each centrality measure with some buffer
        degree_max = maximum(df.degree_centrality) * 1.15
        betweenness_max = maximum(df.betweenness_centrality) * 1.15
        closeness_max = maximum(df.closeness_centrality) * 1.15
        eigenvector_max = maximum(df.eigenvector_centrality) * 1.15
        
        # Store the maximum overall for this network type
        y_limits[nt] = max(degree_max, betweenness_max, closeness_max, eigenvector_max)
    end
    
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
        
        # Set y-axis limits based on the specific network type
        # This helps better visualize the distribution patterns
        y_max = y_limits[nt]
        
        # Create a more precisely aligned boxplot using @df macro
        p = @df boxplot_data boxplot(
            :centrality_type, 
            :value,
            fillcolor=[:blue :red :gray :gold],
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
            xrotation=0,     
            ylims=(0, y_max), # Use network-specific y-axis limit
            ylabel="Centrality Value"
        )
        
        push!(plots, p)
    end
    
    # Combine plots in a single figure with improved margins
    combined_plot = plot(plots..., 
                      layout=plot_layout, 
                      size=(1200, 450),
                      margin=4mm,         
                      left_margin=6mm,    
                      right_margin=2mm,   
                      bottom_margin=8mm,  
                      top_margin=6mm,     
                      xtickfontsize=9,
                      link=:none,        # Don't link y axes to allow different scales
                      title_position=:center)
    
    # Save the figure
    savefig(combined_plot, "figures/centrality_comparison_mdeg_$(mean_degree).pdf")
    
    return combined_plot
end