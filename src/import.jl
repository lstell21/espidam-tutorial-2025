using Graphs, GraphIO.EdgeList

graph = loadgraph("degs/network", "graph_key", EdgeListFormat())

#remove all double edges
for e in edges(graph)
    if has_edge(graph, src(e), dst(e))
        rem_edge!(graph, src(e), dst(e))
    end
end

#convert to undirected
graph = SimpleGraph(graph)

using GraphPlot

gplot(graph)