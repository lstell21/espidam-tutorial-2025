############################## Necessary setup ##############################
begin
    # Navigate to the directory containing this file
    cd(@__DIR__)

    # Create directories if they don't exist
    if !isdir("data")
        mkdir("data")
    end
    if !isdir("figures")
        mkdir("figures")
    end
    if !isdir("output")
        mkdir("output")
    end

    # Import the Pkg module, activate the current directory as the environment for Pkg, instantiate the environment
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()

    # Import the necessary packages
    using Agents, Graphs, Random, Plots, DataFrames, CSV, Statistics, Distributions

    # Read in deg_dist.csv to use for the random_configuration_model
    degrees = CSV.read("deg_dist.csv", DataFrame, header=false)
end
############################## Necessary setup ##############################


################################ Model setup ################################
begin
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
end
################################ Model setup ################################


# Initialization: initialize the model with the chosen parameters: network_type, patient_zero, high_risk, mean_degree
model = initialize(; network_type=:random)

include("src/analyze_graph.jl")

# Analyze the graph, output the results to a CSV file
graph_measures = analyze_graph(model.graph)
CSV.write("data/graph_keyfig_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).csv", graph_measures)

include("src/plotting.jl")

#plot_degree_distribution(graph_measures)
plotdegdist = plot_degree_distribution(graph_measures)
display(plotdegdist)
savefig(plotdegdist, "figures/plotdegdist_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).pdf")

# specify "adata" (agent data to collect) and "mdata" (model data to collect)
adata = [:status]
mdata = [:susceptible_count, :infected_count, :recovered_count]

adf, mdf = run!(model, 100; adata, mdata)

plotdynamics = plot_epidemic_trajectories(mdf)
display(plotdynamics)
savefig(plotdynamics, "figures/plotdynamics_$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hirisk_$(model.high_risk)_hr_frac_$(model.fraction_high_risk)_trans_$(model.trans_prob).pdf")

# for easier comparison, we plot the same scenario for all network types (initializes and runs model)
# plot_single_run(; network_type=:random)
# plot_single_run(; network_type=:smallworld)
# plot_single_run(; network_type=:preferentialattachment)
# plot_single_run(; network_type=:configuration)
# plot_single_run(; network_type=:proportionatemixing)

include("src/run_simulations.jl")

# Run the simulation multiple times with different seeds (100 times), output results to file and plot
run_and_plot(; plot_type=:sfr, network_type=:random)

################################## Part 1 ###################################

# Task 1:

# 1. Choose a graph type random
# 2. Run the model; these steps will be done automatically

#   2.1 Generate network (done with Initialization)
#   2.2 Calculate and visualize degree distribution
#   2.3 Calculate network properties (clustering coefficient, component sizes, diameter)
#   2.4. Calculate centrality measures for nodes
#   2.5 Run SIR model on the network with default parameters

# 3. Import the output files into R or Excel, or any other software you like to work with
# 4. Make figures of temporal dynamics, final size distribution, distribution of durations.
# 5. Start with one graph type to go through all steps.

# Task 2:

# 1. Now repeat the same with graph types smallworld and preferentialattachment
# 2. Compare degree distributions and other graph-based measures
# 3. Compare the dynamics of the outbreaks
# 4. Compare the final sizes
# 5. Make graphs to compare the epidemics on these 3 types of networks
# 6. Investigate whether there is a relationship between any of the graph measures and the final epidemic size (number of infected individuals at the end of the outbreak).

# Task 3:

# 1. Now choose one of the 3 network types
# 2. Vary the mean degree and observe the effect on epidemic dynamics and final state.
# 3. Vary the index case (start the outbreak from an individual with high degree_centrality, high betweenness_centrality or high eigenvector_centrality) and observe how it affects the epidemic.

################################## Part 1 ###################################


################################## Part 2 ###################################
# We will now explore the network type configuration network. We do this with a degree distribution based on data, and with a degree distribution sampled from a negative binomial distribution. You got a file with data from the POLYMOD study, from which we will use the column “cnt_count”, which contains the number of daily contacts for all participants.

# Task 1:

# 1. Make yourself familiar with the data file. Choose a country for which you want to extract a degree distribution.
# 2. You need a list of degrees of length 1000. Sample this from the data for the chosen country and write it on a file deg_dist.csv
# 3. Plot a histogram of the distribution
# 4. Visually fit a negative binomial distribution to the data.

# Task 2:

# 1. You now use the file deg_dist.csv to generate a configuration network and run the model on this network.
# 2. Compare the results to those from the previous practical. How does the structure of this network compare to the other network types?
# 3. Interpretation: in which aspects is the generated network not a realistic model of the real contact network measured by POLYMOD?
# 4. Run also the proportionatemixing network with the parameters for the negative binomial distribution estimated (roughly) from the data. Are the results similar to the results from the configuration network?
# 5. Now vary the mean degree of the negative binomial distribution and observe how this influences network structure and epidemic dynamics. (Optional: vary the dispersion parameter and observe how it influences the dynamics).

# Finally we want to investigate the impact of interventions. Note that an agent is characterized by disease status, time since infection, and a variable called risk. The latter we have not used yet. The model distinguishes between high and low risk, where in the default situation there is no difference between the two. In the current implementation 10% of the population are assumed to be high risk. We can modify the model, such that the population is heterogeneous in the following way:

# Task 3:

# 1. In agent_step.jl the transmission of infection is defined. Increase the transmission probability of high risk individuals by a factor between 1 and 5 (trans_prob = 0.1 in default). Alternatively, add a reduction factor for low risk individuals.
# 2. In model.jl at initialization you can choose whether high risk individuals are random, maxdegree, maxbetweenness, or maxeigenvector. This means that high risk individuals will be place on nodes with highest values of these centrality measures.
# 3. Investigate how the location of the high risk individuals (with a higher transmission probability) influences the epidemic dynamics and final size.
# 4. Think about how to define interventions in this model. An intervention like mask wearing for example, would reduce the transmission probability. An intervention could be targeted to high risk individuals or spread evenly across the population.
# 5. Think about how vaccination could be implemented in the model.
# 6. Think about how contact reduction could be implemented in the model.
# 7. If you have some programming experience, you could try to implement these interventions.
# 8. Study how targeting interventions to high risk individuals compares to an intervention that is applied randomly.

################################## Part 2 ###################################