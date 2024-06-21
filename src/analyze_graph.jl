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

    mean_degree = mean(Graphs.degree(g))
    density = Graphs.density(g)
    clustering_coefficient = global_clustering_coefficient(g)
    assortativity = Graphs.assortativity(g)
    diam = is_connected(g) ? diameter(g) : "Graph is not connected"
    degree_distribution = degree_histogram(g)
    dg_c = degree_centrality(g)
    btwn_c = betweenness_centrality(g)
    clns_c = closeness_centrality(g)
    eig_c = eigenvector_centrality(g)
    cnct_components = connected_components(g)
    comp_lengths = map(length, connected_components(g))
    max_comp_length = maximum(map(length, connected_components(g)))
    max_cliques = maximal_cliques(g)

    # Pad where necessary to ensure equal length
    max_length = maximum([length(dg_c), length(btwn_c), length(clns_c), length(eig_c), length(cnct_components), length(comp_lengths), length(max_comp_length), length(max_cliques)])
    mean_degree = vcat(mean_degree, fill(missing, max_length - length(mean_degree)))
    density = vcat(density, fill(missing, max_length - length(density)))
    clustering_coefficient = vcat(clustering_coefficient, fill(missing, max_length - length(clustering_coefficient)))
    assortativity = vcat(assortativity, fill(missing, max_length - length(assortativity)))
    diam = vcat(diam, fill(missing, max_length - 1))
    degree_distribution = vcat(degree_distribution, fill(missing, max_length - 1))
    dg_c = vcat(dg_c, fill(missing, max_length - length(dg_c)))
    btwn_c = vcat(btwn_c, fill(missing, max_length - length(btwn_c)))
    clns_c = vcat(clns_c, fill(missing, max_length - length(clns_c)))
    eig_c = vcat(eig_c, fill(missing, max_length - length(eig_c)))
    cnct_components = vcat(cnct_components, fill(missing, max_length - length(cnct_components)))
    comp_lengths = vcat(comp_lengths, fill(missing, max_length - length(comp_lengths)))
    max_comp_length = vcat(max_comp_length, fill(missing, max_length - length(max_comp_length)))
    max_cliques = vcat(max_cliques, fill(missing, max_length - length(max_cliques)))

    results = DataFrame(
        mean_degree = mean_degree,
        density = density,
        clustering_coefficient = clustering_coefficient,
        assortativity = assortativity,
        diameter = diam,
        degree_distribution = degree_distribution,
        degree_centrality = dg_c,
        betweenness_centrality = btwn_c,
        closeness_centrality = clns_c,
        eigenvector_centrality = eig_c,
        connected_components = cnct_components,
        component_lengths = comp_lengths,
        max_component_length = max_comp_length,
        maximal_cliques = max_cliques
    )
   return results
end