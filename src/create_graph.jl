"""
    create_graph(; network_type, mean_degree, n_nodes, dispersion = 0.1, β = 0.1, k = 3)

Create a graph based on the specified network type.

# Arguments
- `network_type`: The type of network to create. Can be `:random`, `:smallworld`, `:preferentialattachment`, `:configuration` or `:proportionatemixing`.
- `n_nodes`: The number of nodes in the graph. For :configuration, the number of nodes is fixed to 1000. Default is 1000.
- `mean_degree`: The average degree of the nodes in the graph. For :preferentialattachment, k is used instead, k/2 is the mean degree (for even numbers). Default is 4.
- `dispersion`: The dispersion parameter for the negative binomial distribution, used only when `network_type` is `:proportionatemixing`. Default is 0.1.
- `β`: The rewiring probability, used only when `network_type` is `:smallworld`. Default is 0.1.
- `k`: The number of edges to attach from a new node to existing nodes, used only when `network_type` is `:preferentialattachment`. Default is 3.

# Returns
- `graph`: The created graph.

# Examples
```julia
g = create_graph(; network_type = :random,  mean_degree = 4)
```
"""
function create_graph(; network_type::Symbol, mean_degree::Integer, n_nodes::Integer=1000, dispersion::Float64=0.1, β::Float64=0.1, k::Integer=4)
    if network_type == :random
        graph = Graphs.erdos_renyi(n_nodes, mean_degree / n_nodes)
    elseif network_type == :smallworld
        graph = Graphs.newman_watts_strogatz(n_nodes, mean_degree, β::Float64)
    elseif network_type == :preferentialattachment
        graph = Graphs.barabasi_albert(n_nodes, Int(round(mean_degree / 2)))
    elseif network_type == :configuration
        graph = Graphs.random_configuration_model(1000, degrees[!, 1])
    elseif network_type == :proportionatemixing
        graph = Graphs.expected_degree_graph(rand(NegativeBinomial(mean_degree * dispersion / (1 - dispersion), dispersion), n_nodes))
    end
    return graph
end