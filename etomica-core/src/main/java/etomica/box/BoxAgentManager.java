/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.simulation.Simulation;
import etomica.simulation.SimulationEventManager;
import etomica.simulation.SimulationListener;
import etomica.util.collections.IndexMap;

import java.util.Objects;
import java.util.function.Consumer;
import java.util.function.Function;

/**
 * Acts on behalf of client classes (a BoxAgentSource) to manage
 * agents for each Box in a simulation.  When Box instances are added or removed
 * from the simulation, the agents array (indexed by the Box's index) is updated.
 *
 * @author andrew
 */
public final class BoxAgentManager<E> implements SimulationListener {
    private static final Consumer<Object> NULL_RELEASER = o -> {};

    private final Function<Box, ? extends E> agentSource;
    private final Consumer<? super E> agentReleaser;
    private final SimulationEventManager simEventManager;
    private final IndexMap<E> agents;

    /**
     * Constructs and sets the Simulation containing Boxes to be tracked. This will register
     * this BoxAgentManager as a listener to the Simulation so that it will be notified
     * when Boxes are added or removed. This call also initializes the agents for any Boxes
     * already existing in the Simulation.
     * @param source object that makes the agents.
     * @param sim the simulation using this BoxAgentManager
     */
    public BoxAgentManager(BoxAgentSource<E> source, Simulation sim) {
        this(sim, source::makeAgent, source::releaseAgent);
    }

    public BoxAgentManager(Simulation sim, Function<Box, ? extends E> source) {
        this(sim, source, NULL_RELEASER);

    }

    public BoxAgentManager(Simulation sim, Function<Box, ? extends E> source, Consumer<? super E> releaser) {
        this.agentSource = Objects.requireNonNull(source);
        this.agentReleaser = Objects.requireNonNull(releaser);
        this.simEventManager = Objects.requireNonNull(sim).getEventManager();
        // this will crash if the given sim is in the middle of its constructor
        simEventManager.addListener(this);

        agents = new IndexMap<>(1, 1);
        sim.getBoxes().forEach(box -> agents.put(box.getIndex(), agentSource.apply(box)));
    }

    public IndexMap<E> getAgents() {
        return agents;
    }

    /**
     * @return the agent associated with the given Box
     */
    public E getAgent(Box box) {
        return this.agents.get(box.getIndex());
    }

    /**
     * Associates an Agent with a Box.
     * @param box the Box
     * @param agent the Agent
     */
    public void setAgent(Box box, E agent) {
        this.agents.put(box.getIndex(), agent);
    }

    /**
     * Notifies the BoxAgentManager that it should release all agents and
     * stop listening for events from the simulation.
     */
    public void dispose() {
        // remove ourselves as a listener to the old box
        simEventManager.removeListener(this);
        this.agents.values().forEach(agentReleaser);
    }

    public void simulationBoxAdded(Simulation sim, Box box) {
        this.agents.put(box.getIndex(), agentSource.apply(box));
    }

    public void simulationBoxRemoved(Simulation sim, Box box) {
        E agent = this.agents.remove(box.getIndex());
        agentReleaser.accept(agent);
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
}
