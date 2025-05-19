# ESPIDAM Tutorial 2025

## Installation instructions

1. **Visual Studio Code**
   - Install the latest version of [Visual Studio Code](https://code.visualstudio.com/)

2. **Julia Language**
   - Install the latest version of [Julia](https://julialang.org/downloads/)
   - Alternatively, use [JuliaUp](https://github.com/JuliaLang/juliaup) for better version management

3. **VS Code Extensions**
   - Open VS Code and install the Julia extension
   - You can find it by searching for "Julia" in the Extensions tab

4. **Verify Installation**
   - Open the Julia REPL in VS Code by pressing:
     - Windows: `Alt+J, Alt+O`
     - macOS: `Ctrl+J, Ctrl+O`
     - Linux: `Ctrl+J, Ctrl+O`
   - Alternatively, use the Command Palette (F1 or Ctrl/Cmd+Shift+P) and search for "Julia: Start REPL"
   - If the REPL opens successfully, your setup is working correctly

5. **Run the Tutorial**
   - Clone or download this repository
   - Open the notebook files (Part1.jpynb and Part2.jpynb) in VS Code
   
If you encounter any issues during installation, please check the [Julia documentation](https://docs.julialang.org/) or open an issue in this repository.

## Part 1: Network Types and Epidemic Modeling

### Task 1: Getting Started with Random Networks

1. Start with graph type: `random`
2. Run the notebook; these steps will be done automatically:
   - Generate network (done with initialization)
   - Calculate and visualize degree distribution
   - Calculate network properties (clustering coefficient, component sizes, diameter)
   - Calculate centrality measures for nodes
   - Run SIR model on the network with default parameters
3. Start with one graph type to go through all steps

### Task 2: Comparing Network Types

1. Repeat the same analysis with graph types `smallworld` and `preferential`
2. Compare degree distributions and other graph-based measures
3. Compare the dynamics of the outbreaks
4. Compare the final sizes
5. Investigate whether there is a relationship between any of the graph measures and the final epidemic size (number of infected individuals at the end of the outbreak)

### Task 3: Parameter Variation

1. Choose one of the 3 network types
2. Vary the mean degree and observe the effect on epidemic dynamics and final state
3. Vary the index case (start the outbreak from an individual with high degree_centrality, high betweenness_centrality or high eigenvector_centrality) and observe how it affects the epidemic

---

## Part 2: Configuration Networks and Real-World Data

We will now explore the network type configuration network. We do this with a degree distribution based on data, and with a degree distribution sampled from a negative binomial distribution. You got a file with data from the POLYMOD study, from which we will use the column "cnt_count", which contains the number of daily contacts for all participants.

### Task 1: Data Preparation

1. Make yourself familiar with the data file. Choose a country for which you want to extract a degree distribution
2. You need a list of degrees of length 1000. Sample this from the data for the chosen country and write it on a file `deg_dist.csv`
3. Plot a histogram of the distribution
4. Visually fit a negative binomial distribution to the data

### Task 2: Configuration Networks

1. Use the file `deg_dist.csv` to generate a configuration network and run the model on this network
2. Compare the results to those from the previous practical. How does the structure of this network compare to the other network types?
3. Interpretation: in which aspects is the generated network not a realistic model of the real contact network measured by POLYMOD?
4. Run also the proportionatemixing network with the parameters for the negative binomial distribution estimated (roughly) from the data. Are the results similar to the results from the configuration network?
5. Now vary the mean degree of the negative binomial distribution and observe how this influences network structure and epidemic dynamics. (Optional: vary the dispersion parameter and observe how it influences the dynamics)

### Task 3: Risk Factors and Interventions

1. In `agent_step.jl` the transmission of infection is defined. Increase the transmission probability of high risk individuals by a factor between 1 and 5 (`trans_prob = 0.1` in default). Alternatively, add a reduction factor for low risk individuals.
2. In `model.jl` at initialization you can choose whether high risk individuals are random, maxdegree, maxbetweenness, or maxeigenvector. This means that high risk individuals will be place on nodes with highest values of these centrality measures.
3. Investigate how the location of the high risk individuals (with a higher transmission probability) influences the epidemic dynamics and final size.
4. Think about how to define interventions in this model. An intervention like mask wearing for example, would reduce the transmission probability. An intervention could be targeted to high risk individuals or spread evenly across the population.
5. Think about how vaccination could be implemented in the model.
6. Think about how contact reduction could be implemented in the model.
7. If you have some programming experience, you could try to implement these interventions.
8. Study how targeting interventions to high risk individuals compares to an intervention that is applied randomly.