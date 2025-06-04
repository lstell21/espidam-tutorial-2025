"""
    agent_step!(person::Person, model::ABM)

Update the state of an agent in the ABM model for one time step.

# Arguments
- `person::Person`: The agent whose state needs to be updated.
- `model::ABM`: The ABM model containing the agent.

# Description
- If infected, increment days infected and recover if enough time has passed
- If infected, potentially transmit to nearby susceptible agents with probability
  adjusted by risk level
"""
function agent_step!(person::Person, model::ABM)
    # Handle recovery if infected
    if person.status == :I
        person.days_infected += 1
        if person.days_infected >= model.days_to_recovered
            person.status = :R
            return
        end
        
        # Handle infection of neighbors
        transmission_probability = model.trans_prob
        if person.risk == :low
            transmission_probability *= model.low_risk_factor
        end
        
        for neighbor in nearby_agents(person, model, 1)
            if neighbor.status == :S && rand(abmrng(model)) < transmission_probability
                neighbor.status = :I
            end
        end
    end
end