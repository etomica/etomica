package etomica.atom.iterator;

import etomica.AtomPair;
import etomica.AtomSet;
import etomica.Phase;
import etomica.Species;
import etomica.IteratorDirective.Direction;

/**
 * Iterator that returns a single pair of atoms taken from one or two species
 * according to a target specification. Wraps an instance of Api11 and sets its
 * basis to be the species agents of the species in the phase.
 * 
 * @see Api11
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on May 10, 2005 by kofke
 */
public class ApiSpecies11 extends AtomPairIteratorAdapter implements
        AtomsetIteratorMolecule {

    /**
     * Constructs iterator to iterate molecules from two species. Species may be
     * the same or different instances, and cannot be changed after
     * construction.
     * 
     * @param species
     *            length-2 array of species, which may refer to the same species
     *            instance or to two different species instances
     * @throws IllegalArgumentException
     *             if species.length != 2
     * @throws NullPointerException
     *             if species or either of its elements is null
     */
    public ApiSpecies11(Species[] species) {
        super(new Api11());
        if (species.length != 2)
            throw new IllegalArgumentException(
                    "ApiSpecies11 requires array of two species");
        if (species[0] == null || species[1] == null)
            throw new NullPointerException(
                    "Constructor of ApiInterspecies11 requires two non-null species");
        species0 = species[0];
        species1 = species[1];
        setPhase(null);
    }

    /**
     * Sets the phase from which iterates are taken.
     */
    public void setPhase(Phase phase) {
        if (phase != null) {
            pair.atom0 = phase.getAgent(species0);
            pair.atom1 = phase.getAgent(species1);
            ((Api11) iterator).setBasis(pair);
        }

    }

    /**
     * Performs no action. Implementation of AtomsetIteratorMolecule interface
     */
    public void setDirection(Direction direction) {
        //performs no action
    }

    /**
     * Sets the target atoms that determine the atoms given in the iterate. Atom
     * set must contain at least two atoms. Target must identify exactly two
     * child atoms (one for each species) through which all target atoms are
     * descended. If these conditions are not met, then no iterate will be given
     * on reset.
     * 
     * @throws NullPointerException
     *             if targetAtoms is null
     * @throws IllegalArgumentException
     *             if targetAtoms.count() is < 2
     */
    public void setTarget(AtomSet targetAtoms) {
        ((Api11) iterator).setTarget(targetAtoms);
    }

    private final Species species0, species1;
    private final AtomPair pair = new AtomPair();
}
