"""
    analyze_graph(g)

Analyze the given graph `g` and return a DataFrame containing various graph metrics.

# Arguments
- `g`: The input graph.

# Returns
A DataFrame containing the following graph metrics:
- `"density"`: The density of the graph.
- `"mean_degree"`: The mean degree of the graph.
- `"clustering_coefficient"`: The clustering coefficient of the graph.
- `"assortativity"`: The assortativity of the graph.
- `"connected_components"`: The number of connected components in the graph.
- `"component_lengths"`: The lengths of each connected component in the graph.
- `"max_component_length"`: The length of the largest connected component in the graph.
- `"maximal_cliques"`: The maximal cliques in the graph.
- `"diameter"`: The diameter of the graph if it is connected, otherwise "Graph is not connected".
- `"degree_distribution"`: The degree distribution of the graph.
- `"degree_centrality"`: The degree centrality of the graph.
- `"betweenness_centrality"`: The betweenness centrality of the graph.
- `"closeness_centrality"`: The closeness centrality of the graph.
- `"eigenvector_centrality"`: The eigenvector centrality of the graph.
"""
function analyze_graph(g::AbstractGraph)
    # Calculate network metrics
    mean_degree = mean(Graphs.degree(g))
    density = Graphs.density(g)
    clustering_coefficient = global_clustering_coefficient(g)
    assortativity = Graphs.assortativity(g)
    diam = is_connected(g) ? diameter(g) : "Graph is not connected"
    degree_distribution = degree_histogram(g)
    
    # Calculate centrality measures
    dg_c = degree_centrality(g)
    btwn_c = betweenness_centrality(g)
    clns_c = closeness_centrality(g)
    eig_c = eigenvector_centrality(g)
    
    # Calculate component information
    cnct_components = connected_components(g)
    comp_lengths = map(length, connected_components(g))
    max_comp_length = maximum(map(length, connected_components(g)))
    max_cliques = maximal_cliques(g)
    
    # Create summary dataframe with scalar values
    results_summary = DataFrame(
        metric = [
            "Mean Degree", 
            "Density", 
            "Clustering Coefficient",
            "Assortativity",
            "Diameter",
            "Number of Connected Components",
            "Max Component Length"
        ],
        value = [
            mean_degree, 
            density, 
            clustering_coefficient,
            assortativity,
            diam,
            length(cnct_components),
            max_comp_length
        ]
    )
    
    # Store more detailed node-level metrics in a separate dataframe
    centrality_measures = DataFrame(
        node_id = 1:nv(g),
        degree = Graphs.degree(g),
        degree_centrality = dg_c,
        betweenness_centrality = btwn_c,
        closeness_centrality = clns_c,
        eigenvector_centrality = eig_c
    )
    
    # Create a dictionary to hold all results
    full_results = Dict(
        "summary" => results_summary,
        "centrality" => centrality_measures,
        "degree_distribution" => degree_distribution,
        "connected_components" => cnct_components,
        "component_lengths" => comp_lengths,
        "maximal_cliques" => max_cliques
    )
    
    return full_results
end

"""
    print_graph_analysis(analysis_results)

Display the graph analysis results in a notebook-friendly format.

# Arguments
- `analysis_results`: The dictionary returned by `analyze_graph`.
"""
function print_graph_analysis(analysis_results)
    # Create a result that will display nicely in Jupyter
    summary_df = analysis_results["summary"]
    
    # Get centrality statistics for display
    centrality_summary = describe(analysis_results["centrality"][:, 2:end])
    
    # Get degree distribution info
    deg_dist = analysis_results["degree_distribution"]
    deg_info = DataFrame(
        metric = ["Min Degree", "Max Degree", "Most Common Degree", "Number of Connected Components"],
        value = [
            minimum(keys(deg_dist)), 
            maximum(keys(deg_dist)),
            findmax(collect(values(deg_dist)))[2],
            length(analysis_results["connected_components"])
        ]
    )
    
    # Return a named tuple of dataframes that will display well in Jupyter
    return (
        summary = summary_df,
        centrality = centrality_summary,
        degree_info = deg_info,
        component_sizes = DataFrame(size = analysis_results["component_lengths"])
    )
end