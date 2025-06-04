"""
initialize(; network_type, mean_degree = 4, n_nodes = 1000, dispersion = 0.1, patient_zero = :random, high_risk=:random, fraction_high_risk=0.1, trans_prob = 0.1, days_to_recovered = 14, seed = 42, r̂ = nothing, p̂ = nothing, low_risk_factor = 1.0)

Initialize the model with default parameters.

# Arguments
- `network_type`: The type of network to create. Can be `:random`, `:smallworld`, `:preferentialattachment`, `:configuration` or `:proportionatemixing`.
- `mean_degree`: The mean degree of the network. For :preferentialattachment, k is used as mean_degree/2 (for even numbers). Default is 4.
- `n_nodes`: The number of nodes in the network. For :configuration, the number of nodes is fixed to 1000. Default is 1000.
- `dispersion`: The dispersion parameter for the negative binomial distribution, used only when `network_type` is `:proportionatemixing` and r̂ and p̂ are not provided. Default is 0.1.
- `patient_zero`: The type of patient zero. Can be `:random`(a random agent), `:maxdegree` (the agent with highest degree_centrality), `:maxbetweenness`, and `:maxeigenvector`.
- `high_risk`: The distribution of high and low risk agents. Can be `:random` (randomly distributed), `:maxdegree` (based on degree centrality), `:maxbetweenness` (based on betweenness centrality), and `:maxeigenvector` (based on eigenvector centrality).
- `fraction_high_risk`: The fraction of high risk agents in the network. Default is 0.1.
- `trans_prob`: The transmission probability of the disease. Default is 0.1.
- `days_to_recovered`: The number of days it takes for an agent to recover. Default is 14.
- `seed`: The seed for the random number generator. Default is 42.
- `r̂`: The r parameter for negative binomial distribution, used only when `network_type` is `:proportionate`. If not provided, will be calculated from mean_degree and dispersion.
- `p̂`: The p parameter for negative binomial distribution, used only when `network_type` is `:proportionate`. If not provided, will be calculated from mean_degree and dispersion.
- `low_risk_factor`: Factor to multiply the transmission probability for low risk agents. Default is 1.0.

# Returns
- `model`: The created model.
"""
function initialize(; network_type::Symbol, mean_degree::Integer=4, n_nodes::Integer=1000, dispersion::Float64=0.1, patient_zero::Symbol=:random, high_risk::Symbol=:random, fraction_high_risk::Float64=0.1, trans_prob::Float64=0.1, days_to_recovered::Integer=14, seed=42, r̂=nothing, p̂=nothing, low_risk_factor::Float64=1.0)
    # create a graph space
    graph = create_graph(; network_type, mean_degree, n_nodes, dispersion, r̂, p̂)
    space = GraphSpace(graph)
    # set up properties
    properties = create_properties(graph, network_type, n_nodes, mean_degree, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob, days_to_recovered, low_risk_factor, r̂, p̂)
    # set up RNG
    rng = Xoshiro(seed)
    # create the model
    model = ABM(Person, space; properties, agent_step!, model_step!, rng)
    # add agents, if high_risk is random, add high risk agents randomly
    populate(model, high_risk, fraction_high_risk)
    set_patient_zero!(model, patient_zero)
    return model
end


############################### helper functions ###############################

"""
create_properties(graph, network_type, n_nodes, mean_degree, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob, days_to_recovered, low_risk_factor, r̂, p̂)

Create a dictionary of properties for the simulation.

# Arguments
- `graph`: The graph object representing the network structure.
- `network_type`: The type of network.
- `n_nodes`: The number of nodes in the network.
- `mean_degree`: The mean degree of the network.
- `dispersion`: The dispersion parameter for the network.
- `patient_zero`: The initial infected node.
- `high_risk`: Whether high-risk nodes are present in the network.
- `fraction_high_risk`: The fraction of high-risk nodes in the network.
- `trans_prob`: The transmission probability.
- `days_to_recovered`: The number of days it takes for an infected node to recover.
- `low_risk_factor`: Factor to multiply the transmission probability for low risk agents.
- `r̂`: The r parameter for negative binomial distribution.
- `p̂`: The p parameter for negative binomial distribution.

# Returns
- `properties`: A dictionary containing the properties for the simulation.

"""
function create_properties(graph, network_type, n_nodes, mean_degree, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob, days_to_recovered, low_risk_factor=1.0, r̂=nothing, p̂=nothing)
    properties = Dict(
        :graph => graph,
        :network_type => network_type,
        :n_nodes => n_nodes,
        :mean_degree => mean_degree,
        :dispersion => dispersion,
        :patient_zero => patient_zero,
        :high_risk => high_risk,
        :fraction_high_risk => fraction_high_risk,
        :trans_prob => trans_prob,
        :days_to_recovered => days_to_recovered,
        :low_risk_factor => low_risk_factor,
        :susceptible_count => n_nodes,
        :infected_count => 1,
        :recovered_count => 0)
    
    # Add r̂ and p̂ to properties if provided
    if r̂ !== nothing
        properties[:r̂] = r̂
    end
    if p̂ !== nothing
        properties[:p̂] = p̂
    end
    
    return properties
end

"""
populate(model, high_risk::Symbol, fraction_high_risk::Float64)

Populates the model with agents based on the specified risk distribution.

## Arguments
- `model::ABM`: The model to populate with agents.
- `high_risk::Symbol`: The risk distribution to use. Possible values are `:random`, `:maxdegree`, `:maxbetweenness`, and `:maxeigenvector`.
- `fraction_high_risk::Float64`: The fraction of nodes to assign as high-risk agents. Default is 0.1.

## Details
- If `high_risk` is `:random`, agents are randomly assigned as high-risk or low-risk.
- If `high_risk` is `:maxdegree`, fraction_high_risk agents with the highest degree centrality are assigned as high-risk.
- If `high_risk` is `:maxbetweenness`, fraction_high_risk agents with the highest betweenness centrality are assigned as high-risk.
- If `high_risk` is `:maxeigenvector`, fraction_high_risk agents with the highest eigenvector centrality are assigned as high-risk.
"""
function populate(model::ABM, high_risk::Symbol, fraction_high_risk::Float64=0.1)
    if high_risk == :random
        for _ in 1:Int(fraction_high_risk * model.n_nodes)
            add_agent_single!(model, :S, 0, :high)
        end
        for _ in 1:(model.n_nodes-fraction_high_risk*model.n_nodes)
            add_agent_single!(model, :S, 0, :low)
        end
    elseif high_risk == :maxdegree
        sorted_nodes = sortperm(degree_centrality(model.graph), rev=true)
        selected_positions = sorted_nodes[1:Int(floor(fraction_high_risk * length(sorted_nodes)))]
        for i in 1:model.n_nodes
            if i in selected_positions
                add_agent!(i, model, :S, 0, :high)
            else
                add_agent_single!(model, :S, 0, :low)
            end
        end
    elseif high_risk == :maxbetweenness
        sorted_nodes = sortperm(betweenness_centrality(model.graph), rev=true)
        selected_positions = sorted_nodes[1:Int(floor(fraction_high_risk * length(sorted_nodes)))]
        for i in 1:model.n_nodes
            if i in selected_positions
                add_agent!(i, model, :S, 0, :high)
            else
                add_agent_single!(model, :S, 0, :low)
            end
        end
    elseif high_risk == :maxeigenvector
        sorted_nodes = sortperm(eigenvector_centrality(model.graph), rev=true)
        selected_positions = sorted_nodes[1:Int(floor(fraction_high_risk * length(sorted_nodes)))]
        for i in 1:model.n_nodes
            if i in selected_positions
                add_agent!(i, model, :S, 0, :high)
            else
                add_agent_single!(model, :S, 0, :low)
            end
        end
    end
end

"""
set_patient_zero!(model::ABM, patient_zero::Symbol)

Set the initial infected agent in the model based on the given `patient_zero` strategy.

## Arguments
- `model::ABM`: The agent-based model.
- `patient_zero::Symbol`: The strategy to determine the initial infected agent. Possible values are `:random`, `:maxdegree`, `:maxbetweenness`, and `:maxeigenvector`.

## Details
- If `patient_zero` is `:random`, a random agent in the model will be set as infected.
- If `patient_zero` is `:maxdegree`, the agent with the highest degree centrality in the model's graph will be set as infected.
- If `patient_zero` is `:maxbetweenness`, the agent with the highest betweenness centrality in the model's graph will be set as infected.
- If `patient_zero` is `:maxeigenvector`, the agent with the highest eigenvector centrality in the model's graph will be set as infected.
"""
function set_patient_zero!(model::ABM, patient_zero::Symbol)
    if patient_zero == :random
        random_agent(model).status = :I
    elseif patient_zero == :maxdegree
        model[argmax(degree_centrality(model.graph))].status = :I
    elseif patient_zero == :maxbetweenness
        model[argmax(betweenness_centrality(model.graph))].status = :I
    elseif patient_zero == :maxeigenvector
        model[argmax(eigenvector_centrality(model.graph))].status = :I
    end
end