/*
 * History
 * Created on Dec 30, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.phase.Phase;
import etomica.species.Species;

/**
 * Gives pairs formed from the molecules of a single species in a phase. Species
 * is specified at construction and cannot be changed afterwards.
 */
public class ApiIntraspeciesAA extends AtomPairIteratorAdapter implements
        AtomsetIteratorPhaseDependent {

    /**
     * @param species
     *            species whose molecules will form the pair iterates
     */
    public ApiIntraspeciesAA(Species species) {
        this(new Species[] { species, species });
    }

    /**
     * @param species
     *            array of two non-null elements referencing the same species
     *            instance
     * @throws NullPointerException
     *             if species or one of its elements is null
     * @throws IllegalArgumentException
     *             if species is not a length-2 array or if its elements refer
     *             to different species instances
     */
    public ApiIntraspeciesAA(Species[] species) {
        super(new ApiIntraList());
        if (species.length != 2)
            throw new IllegalArgumentException(
                    "Constructor of ApiIntraspecies1A requires two references to the same species instance");
        if (species == null || species[0] == null || species[1] == null)
            throw new NullPointerException(
                    "Constructor of ApiIntraspecies1A requires two non-null species references to the same instance");
        if (species[0] != species[1])
            throw new IllegalArgumentException(
                    "Constructor of ApiIntraspecies1A requires references to the same species instance");
        this.species = species[0];
        pairIterator = (ApiIntraList) iterator;
    }

    /**
     * Configures iterator to return molecules from the set species in the given
     * phase.
     */
    public void setPhase(Phase phase) {
        if (phase == null) {
            emptyList.clear();
            pairIterator.setList(emptyList);
        } else {
            pairIterator
                    .setList(((AtomTreeNodeGroup) phase.getAgent(species).node).childList);
        }
    }

    private final ApiIntraList pairIterator;
    private final Species species;
    private final AtomList emptyList = new AtomList();
}
