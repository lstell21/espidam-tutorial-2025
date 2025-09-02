# Import the necessary packages
using Pkg
Pkg.activate(".")

using Agents, Graphs, Random, Plots, DataFrames, CSV, CategoricalArrays, Statistics, StatsBase, StatsPlots, Distributions, Measures, Optim
using GraphIO.EdgeList, GraphPlot

graph = loadgraph("degs/network", "graph_key", EdgeListFormat())

#remove all double edges
for e in edges(graph)
    if has_edge(graph, src(e), dst(e))
        rem_edge!(graph, src(e), dst(e))
    end
end

#convert to undirected
graph = SimpleGraph(graph)

gplot(graph)

include("src/create_graph.jl")

# Agent creation: agents of type Person and properties status, days_infected and risk
@agent struct Person(GraphAgent)
    status::Symbol = :S #((S)usceptible, (I)nfected, (R)ecovered)
    days_infected::Int = 0 # number of days since infection
    risk::Symbol = :high # something to differentiate agents (here, high and low risk)
end

include("src/initialize.jl")
include("src/agent_step.jl")

# Model step: keep track of the infection numbers
function model_step!(model::ABM)
    model.susceptible_count = sum([model[i].status == :S for i in 1:nv(model.graph)])
    model.infected_count = sum([model[i].status == :I for i in 1:nv(model.graph)])
    model.recovered_count = sum([model[i].status == :R for i in 1:nv(model.graph)])
end

