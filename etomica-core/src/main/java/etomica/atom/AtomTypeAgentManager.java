/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.simulation.Simulation;
import etomica.simulation.SimulationAtomTypeIndexEvent;
import etomica.simulation.SimulationListener;
import etomica.simulation.SimulationSpeciesEvent;
import etomica.species.ISpecies;
import etomica.util.collections.IndexMap;

import java.util.LinkedHashMap;
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
public final class AtomTypeAgentManager<E> implements SimulationListener {

    private final AgentSource<E> agentSource;
    private final IndexMap<E> agents;
    protected Simulation sim;

    public AtomTypeAgentManager(AgentSource<E> source, Simulation sim) {
        this.agentSource = Objects.requireNonNull(source);
        this.sim = Objects.requireNonNull(sim);
        this.agents = new IndexMap<>();

        sim.getEventManager().addListener(this);
        for (ISpecies species : sim.getSpeciesList()) {
            for (AtomType atomType : species.getAtomTypes()) {
                this.agents.put(atomType.getIndex(), agentSource.makeAgent(atomType));
            }
        }
    }

    public IndexMap<E> getAgents() {
        return agents;
    }

    /**
     * Sets the agent associated with the given atom type to be the given agent.  The AtomType must be from the
     * Simulation.  The AtomType's old agent is not "released".  This should be done manually if needed.
     */
    public void setAgent(AtomType atomType, E newAgent) {
        this.agents.put(atomType.getIndex(), newAgent);
    }

    /**
     * Convenience method to return the agent the given AtomType.  For repeated access to the agents from multiple
     * AtomTypes, it might be faster to use the above getAgents method.
     */
    public E getAgent(AtomType type) {
        return agents.get(type.getIndex());
    }

    /**
     * Releases the agents associated with the given AtomType and its children.
     */
    private void releaseAgents(ISpecies parentType) {
        for (AtomType leafType : parentType.getAtomTypes()) {
            E agent = agents.get(leafType.getIndex());
            if (agent != null) {
                agentSource.releaseAgent(agent, leafType);
                agents.remove(leafType.getIndex());
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
        for (AtomType newType : species.getAtomTypes()) {
            if (!this.agents.containsKey(newType.getIndex())) {
                this.agents.put(newType.getIndex(), this.agentSource.makeAgent(newType));
            }
        }
    }

    public void simulationSpeciesRemoved(SimulationSpeciesEvent e) {
        releaseAgents(e.getSpecies());
    }

    public void simulationAtomTypeIndexChanged(SimulationAtomTypeIndexEvent e) {
        E agent = this.agents.remove(e.getIndex());
        if (agent != null) {
            this.agents.put(e.getAtomType().getIndex(), agent);
        }
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
