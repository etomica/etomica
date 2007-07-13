package etomica.atom.iterator;

import java.io.Serializable;

import etomica.action.AtomsetAction;
import etomica.atom.AtomAddressManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.species.Species;

/**
 * Gives pairs formed from the molecules of a species in a box, taking one
 * molecule the species with all of its other molecules. Species is specified at
 * construction and cannot be changed afterwards.
 */

public class ApiIntraspecies1A extends ApiSequence1A implements
        AtomsetIteratorPDT, Serializable {

    /**
     * @param species
     *            species whose molecules will form the pair iterates
     */
    public ApiIntraspecies1A(Species species) {
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
    public ApiIntraspecies1A(Species[] species) {
        super();
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
    }

    /**
     * Configures iterator to return molecules from the set species in the given
     * box.
     * @throws NullPointerException if the Box is null
     */
    public void setBox(Box box) {
        if (box == null) {
            throw new IllegalArgumentException("You are a bad person.  I didn't even care about the box, but since you passed null, I'm going to quit.");
        }
    }

    /**
     * Sets the target molecule with which all pairs are formed. Molecule is
     * determined from the atom specified by the atomSet (which must have
     * count() == 1), which may be the molecule itself or an atom that is part
     * of it. If the atom is not in the species given at
     * construction, no iterates will be returned.
     * 
     * @throws NullPointerException
     *             if targetAtom is null
     */
    public void setTarget(IAtom newTargetAtom) {
        if (newTargetAtom == null) {
            throw new NullPointerException("target atom must not be null");
        }
        targetAtom = newTargetAtom;
    }
    
    public void allAtoms(AtomsetAction action) {
        reset();
        super.allAtoms(action);
    }
    
    public void reset() {
        identifyTargetMolecule();
        super.reset();
    }

    /**
     * Finds target molecule as indicated by the target atom. Sets target
     * molecule to null if target atom is null, box is null, or atom is not
     * part of either species.
     */
    private void identifyTargetMolecule() {
        if (targetAtom == null) {
            targetMolecule = null;
        } else {
            targetMolecule = targetAtom;
            while (targetMolecule.getType().getDepth() > AtomAddressManager.MOLECULE_DEPTH) {
                targetMolecule = targetMolecule.getParentGroup();
            }
            if (targetMolecule.getType().getSpecies() != species) {
                targetMolecule = null;
            }
        }
        //targetMolecule may be null here
        setAtom(targetMolecule);
    }

    private static final long serialVersionUID = 2L;
    private final Species species;

    private IAtom targetAtom, targetMolecule;

}
