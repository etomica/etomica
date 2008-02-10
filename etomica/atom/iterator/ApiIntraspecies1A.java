package etomica.atom.iterator;

import java.io.Serializable;

import etomica.action.AtomsetAction;
import etomica.atom.AtomToAtomSetSpecies;
import etomica.atom.IAtom;
import etomica.atom.IAtomLeaf;
import etomica.box.Box;
import etomica.species.ISpecies;

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
    public ApiIntraspecies1A(ISpecies species) {
        this(species, new AtomToAtomSetSpecies(species));
    }
    
    protected ApiIntraspecies1A(ISpecies species, AtomToAtomSetSpecies atomToAtomSet) {
        super(new AtomIteratorArrayList(IteratorDirective.Direction.UP, 1, atomToAtomSet, atomToAtomSet),
                new AtomIteratorArrayList(IteratorDirective.Direction.DOWN, 1, atomToAtomSet, atomToAtomSet));
        if (species == null) {
            throw new NullPointerException("Constructor of ApiIntraspecies1A a non-null species");
        }
        this.atomToAtomSet = atomToAtomSet;
        this.species = species;
    }

    /**
     * Configures iterator to return molecules from the set species in the given
     * box.
     * @throws NullPointerException if the Box is null
     */
    public void setBox(Box newBox) {
        if (newBox == null) {
            throw new IllegalArgumentException("You passed a null Box.  Now sit in the corner.");
        }
        atomToAtomSet.setBox(newBox);
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
            if (targetMolecule instanceof IAtomLeaf) {
                targetMolecule = ((IAtomLeaf)targetMolecule).getParentGroup();
            }
            if (targetMolecule.getType().getSpecies() != species) {
                targetMolecule = null;
            }
        }
        //targetMolecule may be null here
        setAtom(targetMolecule);
    }

    private static final long serialVersionUID = 2L;
    private final ISpecies species;

    private IAtom targetAtom, targetMolecule;
    protected Box box;
    protected AtomToAtomSetSpecies atomToAtomSet;
}
