/*
 * History
 * Created on Dec 30, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.Phase;
import etomica.Species;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;

/**
 * Gives pairs formed from the molecules of a single species in a phase.
 * Species is specified at construction and cannot be changed afterwards.
 */
public class ApiIntraspeciesAA extends AtomsetIteratorAdapter implements
        AtomsetIteratorPhaseDependent {

    /**
     * @param species species whose molecules will form the pair iterates
     */
    public ApiIntraspeciesAA(Species species) {
        this(new Species[] {species, species});
    }
    
    /**
     * @param species array of two non-null elements referencing the same species instance
     */
    public ApiIntraspeciesAA(Species[] species) {
        super(new ApiListSimple());
        if(species == null || species.length < 1 || species[0] == null) throw new NullPointerException("Constructor of ApiIntraspeciesAA requires two non-null species references to the same instance");
        if(species[0] != species[1]) throw new IllegalArgumentException("Constructor of ApiIntraspeciesAA requires references to the same species instance");
        this.species = species[0];
        pairIterator = (ApiListSimple)iterator;
    }

    /** 
     * Configures iterator to return molecules from the set species in the given phase.
     */
    public void setPhase(Phase phase) {
        if(phase == null) {
            emptyList.clear();
            pairIterator.setList(emptyList);
        } else {
            pairIterator.setList(((AtomTreeNodeGroup)phase.getAgent(species).node).childList);
        }
    }

    private final ApiListSimple pairIterator;
    private final Species species;
    private final AtomList emptyList = new AtomList();
}
