using Agents
using MetaGraphs
using Graphs
using Random
using Statistics

# Define our agent type with additional properties
@agent struct SocialAgent(GraphAgent)
    energy::Float64
    influence::Float64
    cooperation_level::Float64
end

# Function to create a social network with different connection types
function create_social_network()
    # Create a base graph (e.g., a small world network)
    base_graph = watts_strogatz(10, 4, 0.3)
    
    # Convert to MetaGraph to add metadata
    mg = MetaGraph(base_graph)
    
    # Add different types of connections with metadata
    connection_types = ["friendship", "professional", "family", "casual"]
    
    for edge in edges(mg)
        src_node, dst_node = edge.src, edge.dst
        
        # Randomly assign connection type
        conn_type = rand(connection_types)
        
        # Set different weights based on connection type
        if conn_type == "family"
            weight = 1.0  # Strongest connection
            trust_level = 0.9
        elseif conn_type == "friendship"
            weight = 0.8
            trust_level = 0.7
        elseif conn_type == "professional"
            weight = 0.6
            trust_level = 0.5
        else  # casual
            weight = 0.3
            trust_level = 0.2
        end
        
        # Set edge metadata
        set_prop!(mg, edge, :connection_type, conn_type)
        set_prop!(mg, edge, :weight, weight)
        set_prop!(mg, edge, :trust_level, trust_level)
    end
    
    return mg
end

# Agent interaction function that uses connection metadata
function agent_step!(agent, model)
    # Get all neighboring agents
    neighbor_ids = nearby_ids(agent, model, 1)  # radius = 1 for direct neighbors
    
    for neighbor_id in neighbor_ids
        neighbor = model[neighbor_id]
        
        # Get the edge between this agent and neighbor
        edge_data = get_edge_metadata(agent.pos, neighbor.pos, model)
        
        if edge_data !== nothing
            interact_agents!(agent, neighbor, edge_data, model)
        end
    end
end

# Function to get edge metadata between two positions
function get_edge_metadata(pos1, pos2, model)
    graph = abmspace(model).graph
    
    if has_edge(graph, pos1, pos2)
        edge_props = Dict()
        edge_props[:connection_type] = get_prop(graph, pos1, pos2, :connection_type)
        edge_props[:weight] = get_prop(graph, pos1, pos2, :weight)
        edge_props[:trust_level] = get_prop(graph, pos1, pos2, :trust_level)
        return edge_props
    end
    return nothing
end

# Interaction function based on connection type and weight
function interact_agents!(agent1, agent2, edge_data, model)
    connection_type = edge_data[:connection_type]
    weight = edge_data[:weight]
    trust_level = edge_data[:trust_level]
    
    # Different interactions based on connection type
    if connection_type == "family"
        # Family connections: high energy sharing and cooperation boost
        energy_transfer = 0.1 * weight * agent1.energy
        agent1.energy -= energy_transfer * 0.5  # Less energy loss for family
        agent2.energy += energy_transfer
        
        # Boost cooperation for family connections
        agent2.cooperation_level = min(1.0, agent2.cooperation_level + 0.05 * trust_level)
        
    elseif connection_type == "friendship"
        # Friends: moderate energy sharing and influence exchange
        energy_transfer = 0.05 * weight * agent1.energy
        agent1.energy -= energy_transfer
        agent2.energy += energy_transfer
        
        # Mutual influence based on trust
        influence_exchange = 0.02 * trust_level
        agent1.influence += influence_exchange * agent2.cooperation_level
        agent2.influence += influence_exchange * agent1.cooperation_level
        
    elseif connection_type == "professional"
        # Professional: focus on influence and cooperation
        if agent1.cooperation_level > 0.5
            agent2.influence += 0.03 * weight * trust_level
            agent1.cooperation_level = max(0.0, agent1.cooperation_level - 0.01)
        end
        
    else  # casual connections
        # Casual: minimal interaction, small random effects
        if rand() < 0.3  # 30% chance of interaction
            small_boost = 0.01 * weight
            agent2.energy += small_boost
            agent1.energy -= small_boost * 0.5
        end
    end
    
    # Ensure values stay within bounds
    agent1.energy = max(0.0, min(1.0, agent1.energy))
    agent2.energy = max(0.0, min(1.0, agent2.energy))
    agent1.cooperation_level = max(0.0, min(1.0, agent1.cooperation_level))
    agent2.cooperation_level = max(0.0, min(1.0, agent2.cooperation_level))
end

# Model initialization
function initialize_model(n_agents=10)
    # Create the social network
    social_graph = create_social_network()
    
    # Create the space
    space = GraphSpace(social_graph)
    
    # Create the model
    model = StandardABM(SocialAgent, space; agent_step!)
    
    # Add agents to random positions
    for _ in 1:n_agents
        pos = rand(1:nv(social_graph))
        energy = rand(0.3:0.01:0.8)
        influence = rand(0.1:0.01:0.5)
        cooperation = rand(0.2:0.01:0.7)
        
        add_agent!(pos, model, energy, influence, cooperation)
    end
    
    return model
end

# Function to analyze network connections
function analyze_connections(model)
    graph = abmspace(model).graph
    connection_stats = Dict()
    
    for edge in edges(graph)
        conn_type = get_prop(graph, edge, :connection_type)
        weight = get_prop(graph, edge, :weight)
        trust = get_prop(graph, edge, :trust_level)
        
        if !haskey(connection_stats, conn_type)
            connection_stats[conn_type] = []
        end
        push!(connection_stats[conn_type], (weight=weight, trust=trust))
    end
    
    println("Connection Analysis:")
    for (conn_type, data) in connection_stats
        avg_weight = mean([d.weight for d in data])
        avg_trust = mean([d.trust for d in data])
        count = length(data)
        println("  $conn_type: $count connections, avg weight: $(round(avg_weight, digits=2)), avg trust: $(round(avg_trust, digits=2))")
    end
end

# Example usage
function run_example()
    # Initialize the model
    model = initialize_model(8)
    
    println("Initial state:")
    for agent in allagents(model)
        println("Agent $(agent.id): Energy=$(round(agent.energy, digits=2)), Influence=$(round(agent.influence, digits=2)), Cooperation=$(round(agent.cooperation_level, digits=2))")
    end
    
    # Analyze the network
    analyze_connections(model)
    
    # Run the simulation for a few steps
    println("\nRunning simulation for 5 steps...")
    for _ in 1:5
        step!(model)
    end
    
    println("\nFinal state:")
    for agent in allagents(model)
        println("Agent $(agent.id): Energy=$(round(agent.energy, digits=2)), Influence=$(round(agent.influence, digits=2)), Cooperation=$(round(agent.cooperation_level, digits=2))")
    end
    
    return model
end

# Run the example
model = run_example()