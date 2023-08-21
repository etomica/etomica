/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.box.Box;
import etomica.box.BoxEventListener;
import etomica.box.BoxMoleculeEvent;
import etomica.box.BoxMoleculeIndexEvent;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.util.collections.IndexMap;

import java.util.HashMap;
import java.util.Map;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Stream;

/**
 * MoleculeAgentManager acts on behalf of client classes (an AgentSource) to
 * manage agents in every IMolecule in a box.  When molecules are added or
 * removed from the box, the agents array is updated.
 *
 * @param <E> molecule agent class
 * @author Andrew Schultz
 */
public final class MoleculeAgentManager<E> implements BoxEventListener {
    private static final BiConsumer<Object, IMolecule> NULL_RELEASER = (o, m) -> {
    };

    private final Box box;
    private final SpeciesManager sm;
    private final Function<IMolecule, ? extends E> agentSource;
    private final BiConsumer<? super E, IMolecule> onRelease;
    private final Map<ISpecies, IndexMap<E>> agents;

    public MoleculeAgentManager(SpeciesManager sm, Box box, MoleculeAgentSource<E> source) {
        this(sm, box, source::makeAgent, source::releaseAgent);
    }

    public MoleculeAgentManager(SpeciesManager sm, Box box, Function<IMolecule, ? extends E> agentSource) {
        this(sm, box, agentSource, NULL_RELEASER);
    }

    public MoleculeAgentManager(SpeciesManager sm, Box box, Function<IMolecule, ? extends E> agentSource, BiConsumer<? super E, IMolecule> onRelease) {
        this.sm = sm;
        this.box = box;
        this.agentSource = agentSource;
        this.onRelease = onRelease;
        this.agents = new HashMap<>();

        box.getEventManager().addListener(this);

        // Initialize each species with a map of molecule index to agent
        sm.getSpeciesList().forEach(species -> agents.put(species, new IndexMap<>()));

        // fill in the map with agents from all the molecules
        IMoleculeList molecules = box.getMoleculeList();
        for (int i = 0; i < molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);

            // Make an agent for the molecule and add it to the IndexMap for its species
            E agent = agentSource.apply(molecule);
            if (agent != null) {
                this.agents.get(molecule.getType()).put(molecule.getIndex(), agent);
            }
        }

    }

    public Map<ISpecies, IndexMap<E>> getAgents() {
        return agents;
    }

    public Stream<E> agentStream() {
        // this concatenates the streams for the agent map of each species
        return this.agents.values().stream().flatMap(map -> map.values().stream());
    }

    /**
     * Returns the agent associated with the given IAtom.  The IAtom must be
     * from the Box associated with this instance.
     */
    public E getAgent(IMolecule a) {
        return this.agents.get(a.getType()).get(a.getIndex());
    }

    /**
     * Sets the agent associated with the given atom to be the given agent.
     * The IAtom must be from the Box associated with this instance.  The
     * IAtom's old agent is not released.  This should be done manually if
     * needed.
     */
    public void setAgent(IMolecule a, E newAgent) {
        agents.get(a.getType()).put(a.getIndex(), newAgent);
    }

    /**
     * Convenience method to return the box the Manager is tracking.
     */
    public Box getBox() {
        return box;
    }

    /**
     * Notifies the AtomAgentManager it should disconnect itself as a listener.
     */
    public void dispose() {
        // remove ourselves as a listener to the box
        box.getEventManager().removeListener(this);
        for (ISpecies species : sm.getSpeciesList()) {
            IMoleculeList molecules = box.getMoleculeList(species);
            IndexMap<E> speciesAgents = this.agents.get(species);

            for (int j = 0; j < molecules.size(); j++) {
                IMolecule molecule = molecules.get(j);
                E agent = speciesAgents.remove(molecule.getIndex());
                if (agent != null) {
                    onRelease.accept(agent, molecule);
                }
            }
        }
    }

    public void boxMoleculeAdded(BoxMoleculeEvent e) {
        IMolecule mole = e.getMolecule();
        this.agents.get(mole.getType()).put(mole.getIndex(), this.agentSource.apply(mole));
    }

    public void boxMoleculeRemoved(BoxMoleculeEvent e) {
        IMolecule mole = e.getMolecule();
        IndexMap<E> speciesAgents = this.agents.get(mole.getType());
        E agent = speciesAgents.remove(mole.getIndex());
        if (agent != null) {
            this.onRelease.accept(agent, mole);
        }
    }

    public void boxMoleculeIndexChanged(BoxMoleculeIndexEvent e) {
        IMolecule mole = e.getMolecule();
        // the atom's index changed.  assume it would get the same agent
        int oldIndex = e.getOldIndex();
        IndexMap<E> speciesAgents = this.agents.get(mole.getType());
        E agent = speciesAgents.remove(oldIndex);
        speciesAgents.put(mole.getIndex(), agent);
    }

    /**
     * Interface for an object that wants an agent associated with each
     * IMolecule a Box.
     */
    public interface MoleculeAgentSource<E> {

        /**
         * Returns an agent for the given Atom.
         */
        E makeAgent(IMolecule a);

        /**
         * This informs the agent source that the agent is going away and that
         * the agent source should disconnect the agent from other elements
         */
        void releaseAgent(E agent, IMolecule atom);
    }
}
