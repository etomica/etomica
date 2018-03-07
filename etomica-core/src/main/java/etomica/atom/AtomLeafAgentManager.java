/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.box.*;
import etomica.molecule.IMolecule;
import etomica.util.collections.IndexMap;

import java.util.Objects;
import java.util.function.Function;

/**
 * AtomAgentManager acts on behalf of client classes (an AgentSource) to
 * manage agents for every leaf Atom in a box.  When leaf atoms are added or
 * removed from the box, the agents array (indexed by the atom's global
 * index) is updated.  The client can access and modify the agents via getAgent
 * and setAgent.
 *
 * @author Andrew Schultz
 */
public final class AtomLeafAgentManager<E> extends BoxEventListenerAdapter {

    private final IndexMap<E> agents;
    private final AgentSource<E> agentSource;
    protected final Box box;

    public AtomLeafAgentManager(Function<IAtom, ? extends E> agentSource, Box box) {
        this(new AgentSource<E>() {
            @Override
            public E makeAgent(IAtom a, Box agentBox) {
                return agentSource.apply(a);
            }

            @Override
            public void releaseAgent(E agent, IAtom atom, Box agentBox) {
            }
        }, box);

    }

    public AtomLeafAgentManager(AgentSource<E> source, Box box) {
        Objects.requireNonNull(source);
        Objects.requireNonNull(box);
        this.agentSource = source;
        this.box = box;
        box.getEventManager().addListener(this);

        IAtomList leafList = box.getLeafList();
        agents = new IndexMap<>(leafList.size());
        for(int i = 0; i < leafList.size(); i++) {
            E agent = agentSource.makeAgent(leafList.get(i), box);
            if(agent != null) {
                agents.put(i, agent);
            }
        }
    }

    /**
     * Get the agents map. To iterate over all agents, iterate over the collection
     * returned by agents.values().
     *
     * @return the map from atom leaf index to its agent
     */
    public IndexMap<E> getAgents() {
        return agents;
    }

    /**
     * Returns the agent associated with the given IAtom.  The IAtom must be
     * from the Box associated with this instance.
     */
    public E getAgent(IAtom a) {
        return this.agents.get(a.getLeafIndex());
    }

    public E getAgentUnsafe(int i) {
        return this.agents.getFast(i);
    }

    /**
     * Sets the agent associated with the given atom to be the given agent.
     * The IAtom must be from the Box associated with this instance.  The
     * IAtom's old agent is not released.  This should be done manually if
     * needed.
     */
    public void setAgent(IAtom a, E newAgent) {
        this.agents.put(a.getLeafIndex(), newAgent);
    }

    /**
     * Convenience method to return the box the Manager is tracking.
     */
    public Box getBox() {
        return box;
    }

    /**
     * Notifies the AtomAgentManager it should disconnect itself as a listener and release all of its agents.
     */
    public void dispose() {
        // remove ourselves as a listener to the box
        box.getEventManager().removeListener(this);
        IAtomList leafList = box.getLeafList();
        for(int i = 0; i < leafList.size(); i++) {
            E agent = this.agents.get(i);
            if(agent != null) {
                agentSource.releaseAgent(agent, leafList.get(i), box);
            }
        }
        agents.clear();
    }

    public void boxMoleculeAdded(BoxMoleculeEvent e) {
        IMolecule mole = e.getMolecule();
        // add all leaf atoms below this atom
        IAtomList childList = mole.getChildList();

        for(int i = 0; i < childList.size(); i++) {
            IAtom atom = childList.get(i);
            E agent = this.agentSource.makeAgent(atom, box);
            if(agent != null) {
                this.agents.put(atom.getLeafIndex(), agent);
            }
        }
    }

    public void boxMoleculeRemoved(BoxMoleculeEvent e) {
        IMolecule mole = e.getMolecule();
        // IAtomGroups don't have agents, but nuke all atoms below this atom
        IAtomList childList = mole.getChildList();
        for(int iChild = 0; iChild < childList.size(); iChild++) {
            IAtom childAtom = childList.get(iChild);
            int index = childAtom.getLeafIndex();
            E removed = this.agents.remove(index);
            if(removed != null) {
                agentSource.releaseAgent(removed, childAtom, box);
            }
        }
    }

    public void boxGlobalAtomLeafIndexChanged(BoxIndexEvent e) {
        //TODO: not sure if this is needed anymore
//        // don't use leafList.size() since the SpeciesMaster might be notifying
//        // us that it's about to add leaf atoms
//        int newMaxIndex = e.getIndex();
//        if(agents.length > newMaxIndex + reservoirSize || agents.length < newMaxIndex) {
//            // indices got compacted.  If our array is a lot bigger than it
//            // needs to be, shrink it.
//            // ... or we've been notified that atoms are about to get added to the
//            // system.  Make room for them
//            agents = Arrays.copyOf(agents, newMaxIndex + 1 + reservoirSize);
//        }
    }

    public void boxAtomLeafIndexChanged(BoxAtomIndexEvent e) {
        IAtom a = e.getAtom();
        // the atom's index changed.  assume it would get the same agent
        int oldIndex = e.getIndex();
        E agent = this.agents.get(oldIndex);
        if(agent != null) {
            this.agents.put(a.getLeafIndex(), agent);
            this.agents.remove(oldIndex);
        }
    }

    /**
     * Interface for an object that wants an agent associated with each Atom in
     * a Box.
     */
    public interface AgentSource<E> {

        /**
         * Returns an agent for the given Atom.
         * The agentBox is provided for convenience for agent sources that
         * handle multiple boxes.
         *
         * @param agentBox TODO
         */
        E makeAgent(IAtom a, Box agentBox);

        /**
         * This informs the agent source that the agent is going away and that
         * the agent source should disconnect the agent from other elements.
         * The agentBox is provided for convenience for agent sources that
         * handle multiple boxes.
         */
        void releaseAgent(E agent, IAtom atom, Box agentBox);
    }
}
