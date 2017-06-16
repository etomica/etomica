/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.simulation.*;
import etomica.species.ISpecies;
import etomica.util.Arrays;

import java.io.Serializable;
import java.lang.reflect.Array;

/**
 * AtomTypeAgentManager acts on behalf of client classes (an AgentSource) to 
 * manage agents in every AtomType in a box.  When species are added or 
 * removed from the simulation, the agents array (indexed by the AtomType's 
 * global index) is updated.  The client should call getAgents() at any point 
 * where an atom might have have been added to the system because the old array
 * would be stale at that point.
 * @author andrew
 */
public class SpeciesAgentManager implements SimulationListener, java.io.Serializable {

    private static final long serialVersionUID = 1L;
    private final AgentSource agentSource;
    protected Object[] agents;
    protected Simulation sim;
    
    public SpeciesAgentManager(AgentSource source) {
        agentSource = source;
    }
    
    public SpeciesAgentManager(AgentSource source, Simulation sim) {
        agentSource = source;
        init(sim);
    }
    
    /**
     * Returns an iterator that returns each non-null agent
     */
    public AgentIterator makeIterator() {
        return new AgentIterator(this);
    }
    
    /**
     * Sets the agent associated with the given atom type to be the given
     * agent.  The AtomType must be from the Simulation.  The AtomType's old
     * agent is not "released".  This should be done manually if needed.
     */
    public void setAgent(ISpecies atomType, Object newAgent) {
        agents[atomType.getIndex()] = newAgent;
    }
    
    /**
     * Convenience method to return the agent the given AtomType.  For repeated
     * access to the agents from multiple AtomTypes, it might be faster to use
     * the above getAgents method.
     */
    public Object getAgent(ISpecies type) {
        return agents[type.getIndex()];
    }
    
    /**
     * Releases the agents associated with the given AtomType and its children.
     */
    private void releaseAgents(ISpecies parentType) {
        Object agent = agents[parentType.getIndex()];
        if (agent != null) {
            agentSource.releaseAgent(agent, parentType);
        }
        agents[parentType.getIndex()] = null;
    }
    
    private void makeAllAgents() {
        for (int i=0; i<sim.getSpeciesCount(); i++) {
            addAgent(sim.getSpecies(i));
        }
    }
    
    /**
     * Returns the max index of all the children of the given AtomType
     */
    private int getGlobalMaxIndex() {
        int max = 0;
        for (int i=0; i<sim.getSpeciesCount(); i++) {
            if (sim.getSpecies(i).getIndex() > max) {
                max = sim.getSpecies(i).getIndex();
            }
        }
        return max;
    }
    
    /**
     * Unregisters this class as a listener for AtomType-related events and
     * releases its agents.
     */
    public void dispose() {
        // remove ourselves as a listener to the old box
        sim.getEventManager().removeListener(this);
        for (int i=0; i<sim.getSpeciesCount(); i++) {
            releaseAgents(sim.getSpecies(i));
        }
        agents = null;
    }
    
    /**
     * Sets the SpeciesRoot for which this AtomAgentManager will manage
     * AtomType agents.
     */
    public void init(Simulation newSim) {
        sim = newSim;
        sim.getEventManager().addListener(this);

        int numTypes = getGlobalMaxIndex()+1;

        agents = (Object[])Array.newInstance(agentSource.getSpeciesAgentClass(), numTypes);
        // fill in the array with agents from all the atoms
        makeAllAgents();
    }

    public void simulationSpeciesAdded(SimulationSpeciesEvent e) {
        ISpecies species = e.getSpecies();
        agents = Arrays.resizeArray(agents, species.getIndex()+1);
        addAgent(species);
    }

    public void simulationSpeciesRemoved(SimulationSpeciesEvent e) {
        releaseAgents(e.getSpecies());
    }

    public void simulationSpeciesIndexChanged(SimulationSpeciesIndexEvent e) {
        ISpecies species = e.getSpecies();
        int oldIndex = e.getIndex();
        int newIndex = species.getIndex();
        if (newIndex >= agents.length) {
            agents = Arrays.resizeArray(agents, newIndex+1);
        }
        agents[newIndex] = agents[oldIndex];
        agents[oldIndex] = null;
    }

    public void simulationAtomTypeMaxIndexChanged(SimulationIndexEvent e) {
        int maxIndex = e.getIndex();
        agents = Arrays.resizeArray(agents, maxIndex+1);
    }

    public void simulationBoxAdded(SimulationBoxEvent e) {
    }

    public void simulationBoxRemoved(SimulationBoxEvent e) {
    }

    public void simulationSpeciesMaxIndexChanged(SimulationIndexEvent e) {
    }

    public void simulationAtomTypeIndexChanged(SimulationAtomTypeIndexEvent e) {
    }

    protected void addAgent(ISpecies type) {
        agents[type.getIndex()] = agentSource.makeAgent(type);
    }
    /**
     * Interface for an object that wants an agent associated with each
     * AtomType in a Simulation.
     */
    public interface AgentSource {
        /**
         * Returns the Class of the agent.  This is used to create an array of
         * the appropriate Class.
         */
        Class getSpeciesAgentClass();

        /**
         * Returns an agent for the given AtomType.
         */
        Object makeAgent(ISpecies type);

        /**
         * This informs the agent source that the agent is going away and that
         * the agent source should disconnect the agent from other elements.
         */
        void releaseAgent(Object agent, ISpecies type);
    }

    /**
     * Iterator that loops over the agents, skipping null elements
     */
    public static class AgentIterator implements Serializable {
        private static final long serialVersionUID = 1L;
        private final SpeciesAgentManager agentManager;
        private int cursor;
        private Object[] agents;
        
        protected AgentIterator(SpeciesAgentManager agentManager) {
            this.agentManager = agentManager;
        }

        public void reset() {
            cursor = 0;
            agents = agentManager.agents;
        }

        public boolean hasNext() {
            while (cursor < agents.length) {
                if (agents[cursor] != null) {
                    return true;
                }
                cursor++;
            }
            return false;
        }

        public Object next() {
            cursor++;
            while (cursor-1 < agents.length) {
                if (agents[cursor-1] != null) {
                    return agents[cursor-1];
                }
                cursor++;
            }
            return null;
        }
    }
}
