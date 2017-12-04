/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.simulation.*;
import etomica.util.Arrays;

import java.lang.reflect.Array;

/**
 * Acts on behalf of client classes (a BoxAgentSource) to manage
 * agents for each Box in a simulation.  When Box instances are added or removed
 * from the simulation, the agents array (indexed by the Box's index) is updated.
 *
 * @author andrew
 */
public class BoxAgentManager<E> implements SimulationListener {

    protected final BoxAgentSource<E> agentSource;
    protected SimulationEventManager simEventManager;
    protected E[] agents;

    /**
     * Constructs and sets the Simulation containing Boxes to be tracked. This will register
     * this BoxAgentManager as a listener to the Simulation so that it will be notified
     * when Boxes are added or removed. This call also initializes the agents for any Boxes
     * already existing in the Simulation.
     * @param source object that makes the agents.
     * @param boxAgentClass class of the agent returned by the boxAgentSource
     * @param sim the simulation using this BoxAgentManager
     */
    public BoxAgentManager(BoxAgentSource<E> source, Class boxAgentClass, Simulation sim) {
        agentSource = source;
        simEventManager = sim.getEventManager();
        // this will crash if the given sim is in the middle of its constructor
        simEventManager.addListener(this);

        int boxCount = sim.getBoxCount();
        agents = (E[]) Array.newInstance(boxAgentClass, boxCount);
        if (agentSource == null) {
            return;
        }
        for (int i = 0; i < boxCount; i++) {
            addAgent(sim.getBox(i));
        }
    }

    /**
     * @return the agent associated with the given Box
     */
    public E getAgent(Box box) {
        if (box.getIndex() >= agents.length) return null;
        return agents[box.getIndex()];
    }

    /**
     * Associates an Agent with a Box.
     * @param box the Box
     * @param agent the Agent
     */
    public void setAgent(Box box, E agent) {
        int idx = box.getIndex();
        if (idx >= agents.length) {
            // no room in the array.  reallocate the array with an extra cushion.
            agents = (E[]) Arrays.resizeArray(agents, idx + 1);
        }
        agents[box.getIndex()] = agent;
    }

    /**
     * @return an iterator that returns each non-null agent
     */
    public AgentIterator<E> makeIterator() {
        return new AgentIterator<E>(this);
    }

    /**
     * Notifies the BoxAgentManager that it should release all agents and
     * stop listening for events from the simulation.
     */
    public void dispose() {
        if (agents == null) return;
        // remove ourselves as a listener to the old box
        simEventManager.removeListener(this);
        for (int i = 0; i < agents.length; i++) {
            if (agents[i] != null) {
                agentSource.releaseAgent(agents[i]);
            }
        }
        agents = null;
    }

    public void simulationBoxAdded(SimulationBoxEvent e) {
        addAgent(e.getBox());
    }

    public void simulationBoxRemoved(SimulationBoxEvent e) {
        Box box = e.getBox();
        // The given Box got removed.  The remaining boxes got shifted
        // down.
        int index = box.getIndex();
        if (agentSource != null) {
            agentSource.releaseAgent(agents[index]);
        }
        for (int i = index; i < agents.length - 1; i++) {
            agents[i] = agents[i + 1];
        }
        agents = (E[]) Arrays.resizeArray(agents, agents.length - 1);
    }

    public void simulationSpeciesAdded(SimulationSpeciesEvent e) {
    }

    public void simulationSpeciesRemoved(SimulationSpeciesEvent e) {
    }

    public void simulationSpeciesIndexChanged(SimulationSpeciesIndexEvent e) {
    }

    public void simulationSpeciesMaxIndexChanged(SimulationIndexEvent e) {
    }

    public void simulationAtomTypeIndexChanged(SimulationAtomTypeIndexEvent e) {
    }

    public void simulationAtomTypeMaxIndexChanged(SimulationIndexEvent e) {
    }

    protected void addAgent(Box box) {
        agents = (E[]) Arrays.resizeArray(agents, box.getIndex() + 1);
        if (agentSource != null) {
            agents[box.getIndex()] = agentSource.makeAgent(box);
        }
    }
    /**
     * Interface for an object that makes an agent to be associated with each box
     * upon construction.
     */
    public interface BoxAgentSource<E> {
        /**
         * @param box the Box for this agent
         * @return the object that source wants to have with the given box
         */
        E makeAgent(Box box);

        /**
         * Disconnects the given agent from its box and performs any actions needed to clean up
         * @param agent the agent being released
         */
        void releaseAgent(E agent);
    }

    /**
     * Iterator that loops over the agents, skipping null elements
     */
    public static class AgentIterator<E> {
        private final BoxAgentManager<E> agentManager;
        private int cursor;
        private E[] agents;

        /**
         * @param agentManager BoxAgentManager for the BoxAgent being iterated
         */
        protected AgentIterator(BoxAgentManager<E> agentManager) {
            this.agentManager = agentManager;
        }

        /**
         * Positions iterator to begin the iteration
         */
        public void reset() {
            cursor = 0;
            agents = agentManager.agents;
        }

        /**
         * @return true if the iterator has another agent
         */
        public boolean hasNext() {
            while (cursor < agents.length) {
                if (agents[cursor] != null) {
                    return true;
                }
                cursor++;
            }
            return false;
        }

        /**
         * @return next agent from the iterator
         */
        public E next() {
            cursor++;
            while (cursor - 1 < agents.length) {
                if (agents[cursor - 1] != null) {
                    return agents[cursor - 1];
                }
                cursor++;
            }
            return null;
        }
    }
}
