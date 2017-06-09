/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.simulation.*;
import etomica.util.Arrays;

import java.lang.reflect.Array;

/**
 * BoxAgentManager acts on behalf of client classes (a BoxAgentSource) to manage
 * agents for each Box in a simulation.  When Box instances are added or removed
 * from the simulation, the agents array (indexed by the box's index) is updated.
 * The client should call getAgents() at any point where a Box might have have been
 * added to (or removed from) the system because the old array would be stale at that
 * point.
 *
 * @author andrew
 */
public class BoxAgentManager<E> implements SimulationListener {

    protected final BoxAgentSource<E> agentSource;
    protected final Class boxAgentClass;
    protected SimulationEventManager simEventManager;
    protected E[] agents;

    public BoxAgentManager(BoxAgentSource<E> source, Class boxAgentClass) {
        agentSource = source;
        this.boxAgentClass = boxAgentClass;
        if (source == null) {
            agents = (E[]) Array.newInstance(boxAgentClass, 0);
        }
    }

    public BoxAgentManager(BoxAgentSource<E> source, Class boxAgentClass, Simulation sim) {
        agentSource = source;
        this.boxAgentClass = boxAgentClass;
        setSimulation(sim);
    }

    /**
     * Returns the agent associated with the given box
     */
    public E getAgent(Box box) {
        if (box.getIndex() >= agents.length) return null;
        return agents[box.getIndex()];
    }

    public void setAgent(Box box, E agent) {
        int idx = box.getIndex();
        if (idx >= agents.length) {
            // no room in the array.  reallocate the array with an extra cushion.
            agents = (E[]) Arrays.resizeArray(agents, idx + 1);
        }
        agents[box.getIndex()] = agent;
    }

    /**
     * Returns an iterator that returns each non-null agent
     */
    public AgentIterator<E> makeIterator() {
        return new AgentIterator<E>(this);
    }

    /**
     * Sets the Simulation containing Boxs to be tracked.  This method should
     * not be called if setSimulationEventManager is called.
     */
    public void setSimulation(Simulation sim) {
        simEventManager = sim.getEventManager();
        // this will crash if the given sim is in the middle of its constructor
        simEventManager.addListener(this);

        // hope the class returns an actual class with a null Atom and use it to construct
        // the array
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
     * Notifies the BoxAgentManager that it should release all agents and
     * stop listening for events from the simulation.
     */
    public void dispose() {
        // remove ourselves as a listener to the old box
        if (simEventManager != null) {
            simEventManager.removeListener(this);
        }
        if (agentSource != null) {
            for (int i = 0; i < agents.length; i++) {
                if (agents[i] != null) {
                    agentSource.releaseAgent(agents[i]);
                }
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
     * Interface for an object that makes an agent to be placed in each atom
     * upon construction.  AgentSource objects register with the AtomFactory
     * the produces the atom.
     */
    public interface BoxAgentSource<E> {
        E makeAgent(Box box);

        //allow any agent to be disconnected from other elements
        void releaseAgent(E agent);
    }

    /**
     * Iterator that loops over the agents, skipping null elements
     */
    public static class AgentIterator<E> {
        private final BoxAgentManager<E> agentManager;
        private int cursor;
        private E[] agents;

        protected AgentIterator(BoxAgentManager<E> agentManager) {
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
