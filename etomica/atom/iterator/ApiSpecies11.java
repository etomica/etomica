package etomica.atom.iterator;

import etomica.AtomPair;
import etomica.AtomSet;
import etomica.Phase;
import etomica.Species;
import etomica.IteratorDirective.Direction;


/**
 * Iterator that returns a single pair of atoms taken from
 * one or two species according to a target specification.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 10, 2005 by kofke
 */
public class ApiSpecies11 extends AtomPairIteratorAdapter implements
        AtomsetIteratorMolecule {

    /**
     * @param species length-2 array of species, which may refer to the same species 
     * instance or to two different species instances
     */
    public ApiSpecies11(Species[] species) {
        super(new Api11());
        if(species.length != 2) throw new IllegalArgumentException("ApiSpecies11 requires array of two species");
        if(species[0] == null || species[1] == null) throw new NullPointerException("Constructor of ApiInterspecies11 requires two non-null species");
        species0 = species[0];
        species1 = species[1];
        setPhase(null);
    }

    /* (non-Javadoc)
     * @see etomica.atom.iterator.AtomsetIteratorPhaseDependent#setPhase(etomica.Phase)
     */
    public void setPhase(Phase phase) {
        if(phase != null) {
            pair.atom0 = phase.getAgent(species0);
            pair.atom1 = phase.getAgent(species1);
            ((Api11)iterator).setBasis(pair);
        }

    }

    /**
     * Performs no action.
     */
    public void setDirection(Direction direction) {
        //performs no action
    }

    /**
     * 
     * @throws NullPointerException
     *          if targetAtoms is null
     * @throws IllegalArgumentException
     *          if targetAtoms.count() is < 2
     */
    public void setTarget(AtomSet targetAtoms) {
        ((Api11)iterator).setTarget(targetAtoms);
    }

    private final Species species0, species1;
    private final AtomPair pair = new AtomPair();
}
