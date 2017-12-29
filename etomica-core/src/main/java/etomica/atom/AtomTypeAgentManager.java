/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.simulation.*;
import etomica.species.ISpecies;

import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

/**
 * AtomTypeAgentManager acts on behalf of client classes (an AgentSource) to manage agents for every AtomType in a box.
 * When species are added or removed from the simulation, the agents map is updated.
 * <p>
 * This class should only be used instead of a {@link java.util.Map} when the SimulationListener functionality is
 * required. That is, when the map needs to be kept up-to-date with every AtomType and an object provided by the
 * AgentSource.
 *
 * @author andrew
 */
public class AtomTypeAgentManager<E> implements SimulationListener {

    private final AgentSource<E> agentSource;
    private final Map<AtomType, E> agents;
    protected Simulation sim;

    public AtomTypeAgentManager(AgentSource<E> source, Simulation sim) {
        this.agentSource = Objects.requireNonNull(source);
        this.sim = Objects.requireNonNull(sim);
        this.agents = new HashMap<>();

        sim.getEventManager().addListener(this);
        for (ISpecies species : sim.getSpeciesList()) {
            for (int i = 0; i < species.getAtomTypeCount(); i++) {
                AtomType atomType = species.getAtomType(i);
                this.agents.computeIfAbsent(atomType, agentSource::makeAgent);
            }
        }
    }

    public Map<AtomType, E> getAgents() {
        return agents;
    }

    /**
     * Sets the agent associated with the given atom type to be the given agent.  The AtomType must be from the
     * Simulation.  The AtomType's old agent is not "released".  This should be done manually if needed.
     */
    public void setAgent(AtomType atomType, E newAgent) {
        this.agents.put(atomType, newAgent);
    }

    /**
     * Convenience method to return the agent the given AtomType.  For repeated access to the agents from multiple
     * AtomTypes, it might be faster to use the above getAgents method.
     */
    public E getAgent(AtomType type) {
        return agents.get(type);
    }

    /**
     * Releases the agents associated with the given AtomType and its children.
     */
    private void releaseAgents(ISpecies parentType) {
        for (int i = 0; i < parentType.getAtomTypeCount(); i++) {
            AtomType leafType = parentType.getAtomType(i);
            E agent = agents.get(leafType);
            if (agent != null) {
                agentSource.releaseAgent(agent, leafType);
                agents.remove(leafType);
            }
        }
    }

    /**
     * Unregisters this class as a listener for AtomType-related events and releases its agents.
     */
    public void dispose() {
        sim.getEventManager().removeListener(this);
        for (int i = 0; i < sim.getSpeciesCount(); i++) {
            releaseAgents(sim.getSpecies(i));
        }
    }

    public void simulationSpeciesAdded(SimulationSpeciesEvent e) {
        ISpecies species = e.getSpecies();
        for (int i = 0; i < species.getAtomTypeCount(); i++) {
            AtomType newType = species.getAtomType(i);
            this.agents.computeIfAbsent(newType, this.agentSource::makeAgent);
        }
    }

    public void simulationSpeciesRemoved(SimulationSpeciesEvent e) {
        releaseAgents(e.getSpecies());
    }

    /**
     * Interface for an object that wants an agent associated with each AtomType in a Simulation.
     */
    public interface AgentSource<E> {

        /**
         * Returns an agent for the given AtomType.
         */
        E makeAgent(AtomType type);

        /**
         * This informs the agent source that the agent is going away and that the agent source should disconnect the
         * agent from other elements.
         */
        void releaseAgent(E agent, AtomType type);
    }
}
