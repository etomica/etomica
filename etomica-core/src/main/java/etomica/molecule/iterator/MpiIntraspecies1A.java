/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule.iterator;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeToMoleculeListSpecies;
import etomica.potential.IteratorDirective;
import etomica.species.ISpecies;

import java.io.Serializable;

/**
 * Gives pairs formed from the molecules of a species in a box, taking one
 * molecule the species with all of its other molecules. Species is specified at
 * construction and cannot be changed afterwards.
 */
public class MpiIntraspecies1A implements MoleculesetIteratorPDT, Serializable {

    /**
     * @param species
     *            species whose molecules will form the pair iterates
     */
    public MpiIntraspecies1A(ISpecies species) {
        this(species, new MoleculeToMoleculeListSpecies(species));
    }
    
    protected MpiIntraspecies1A(ISpecies species, MoleculeToMoleculeListSpecies atomToAtomSet) {
        if (species == null) {
            throw new NullPointerException("Constructor of ApiIntraspecies1A requires a non-null species");
        }
        
        aiOuter = new MoleculeIteratorSinglet();
        apiUp = new MpiInnerVariable(aiOuter, new MoleculeIteratorArrayList(IteratorDirective.Direction.UP, 1, atomToAtomSet, atomToAtomSet), false);
        apiDown = new MpiInnerVariable(aiOuter, new MoleculeIteratorArrayList(IteratorDirective.Direction.DOWN, 1, atomToAtomSet, atomToAtomSet), true);

        this.atomToAtomSet = atomToAtomSet;
        this.species = species;
    }

    public int nBody() {
        return 2;
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
    public void setTarget(IMolecule newTargetAtom) {
        if (newTargetAtom == null) {
            throw new NullPointerException("target atom must not be null");
        }
        targetMolecule = newTargetAtom;
    }

    public void reset() {
        // this will throw an exception if you don't call setTarget
        if (targetMolecule.getType() != species) {
            // wrong species.  set ourselves to return nothing
            unset();
            return;
        }
        aiOuter.setMolecule(targetMolecule);

        //upList if specified by direction
        upListNow = (direction != IteratorDirective.Direction.DOWN);
        //dnList only if not explicitly directed up
        doGoDown = (direction != IteratorDirective.Direction.UP);
        
        if (upListNow) {
            apiUp.reset();
        }
        else {
            apiDown.reset();
        }
    }

    public void unset() {
        apiDown.unset();
        upListNow = false;
    }
    
    /**
     * Specifies the direction, which applies only if iterating pairs
     * with a target atom; otherwise, if all pairs from group are indicated,
     * direction is ignored.
     */
    public void setDirection(IteratorDirective.Direction direction) {
        this.direction = direction;
    }
    
    /**
     * Returns the number of atom pairs the iterator will return if
     * reset and iterated in its present state.
     */
    public int size() {
        if (aiOuter.getMolecule() == null) {
            return 0;
        }
        int count = 0;
        if (direction != IteratorDirective.Direction.DOWN) {
            count += apiUp.size();
        }
        if (direction != IteratorDirective.Direction.UP) {
            count += apiDown.size();
        }
        return count;
    }
    
    public IMoleculeList next() {
        if (upListNow) {
            IMoleculeList next = apiUp.next();
            if(next != null || !doGoDown) {
                return next;
            }
            upListNow = false;
            apiDown.reset();
        }
        return apiDown.next();
    }
    
    private static final long serialVersionUID = 3L;
    protected final MoleculeIteratorSinglet aiOuter;
    protected IteratorDirective.Direction direction;
    protected final MpiInnerVariable apiUp, apiDown;
    protected boolean doGoDown, upListNow;

    private final ISpecies species;

    private IMolecule targetMolecule;
    protected Box box;
    protected MoleculeToMoleculeListSpecies atomToAtomSet;
}
