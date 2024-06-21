"""
        agent_step!(person::Person, model::ABM)

Update the state of an agent in the ABM model for one time step.

# Arguments
- `person::Person`: The agent whose state needs to be updated.
- `model::ABM`: The ABM model containing the agent.

# Description
- If the agent is infected (`person.status == :I`), increment the number of days infected (`person.days_infected`) by 1.
    - If the number of days infected is greater than or equal to the number of days required for recovery (`model.days_to_recovered`), change the agent's status to recovered (`person.status = :R`).
- If the agent's risk is high (`person.risk == :high`) and the agent is infected (`person.status == :I`), infect nearby susceptible neighbors with a probability of transmission (`model.trans_prob`).
- If the agent's risk is low (`person.risk == :low`) and the agent is infected (`person.status == :I`), infect nearby susceptible neighbors with a probability of transmission (`model.trans_prob`).

# Example
```julia
agent_step!(person, model)
```
"""
function agent_step!(person::Person, model::ABM)
    # recover if infected for 14 days
    if person.status == :I
        person.days_infected += 1
        if person.days_infected >= model.days_to_recovered
            person.status = :R
        end
    end
    # infect neighbors with probability trans_prob
    if person.risk == :high
        if person.status == :I
            for neighbor in nearby_agents(person, model, 1)
                if neighbor.status == :S
                    if rand(abmrng(model)) < model.trans_prob
                        neighbor.status = :I
                    end
                end
            end
        end
    elseif person.risk == :low
        #add different behavior for low risk agents
        if person.status == :I
            for neighbor in nearby_agents(person, model, 1)
                if neighbor.status == :S
                    if rand(abmrng(model)) < model.trans_prob
                        neighbor.status = :I
                    end
                end
            end
        end
    end
end