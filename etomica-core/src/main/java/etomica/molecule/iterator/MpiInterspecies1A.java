/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule.iterator;

import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IteratorDirective;
import etomica.potential.IteratorDirective.Direction;

import java.io.Serializable;

/**
 * Gives pairs formed from the molecules of two different species in a box,
 * taking one molecule of one species with all molecules of the other. Species
 * are specified at construction and cannot be changed afterwards. The
 * 1-molecule species is identified via the setTarget method, and may be changed
 * from one use of the iterator to the next.
 * <p>
 * Atoms in iterate pairs are ordered such that the first atom is of the first
 * species given in the constructor, and the second atom is of the second
 * species given in the constructor, regardless of the specification of the
 * target.
 * <p>
 * Direction may be specified, and is interpreted as follows:
 * <ul>
 * <li>For null direction, no restriction on iterates.
 * <li>For UP direction, iterates are given if index of target-atom species is
 * less than index of other species.
 * <li>For DOWN direction, iterates are given if index of target-atom species
 * is greater than index of other species.
 * </ul>
 */

public class MpiInterspecies1A implements MoleculesetIteratorPDT,
        Serializable {

    /**
     * Sorts given array of species according to species index, then constructs iterator 
     * such that atom0 of the pair iterates is in (sorted) species[0], and atom1 is in 
     * (sorted) species[1], regardless of which is specified via setTarget.  Thus the
     * species index of atom0 is less than that of atom1, for all iterates.
     * 
     * @param species
     *            array of two different, non-null species
     * 
     * @throws IllegalArgumentException
     *             is species array is not of length 2 or if species in array
     *             refer to the same instance
     * @throws NullPointerException
     *             if species array is null or if either species in array is
     *             null
     */
    public MpiInterspecies1A(ISpecies[] species) {
        super();
        if (species.length != 2) {
            throw new IllegalArgumentException(
                    "Constructor of ApiInterspecies1A requires an array of two species");
        }
        if (species[0] == null || species[1] == null) {
            throw new NullPointerException(
                    "Constructor of ApiInterspeciesAA requires two non-null species");
        }
        if (species[0] == species[1]) {
            throw new IllegalArgumentException(
                    "Constructor of ApiInterspeciesAA requires two different species");
        }
        aiOuter = new MoleculeIteratorSinglet();
        aiInner = new MoleculeIteratorArrayListSimple();
        apiUp = new MpiInnerFixed(aiOuter, aiInner, false);
        apiDown = new MpiInnerFixed(aiOuter, aiInner, true);
        iterator = apiUp;

        // we need to sort these.  we'll do that once we have the box
        species0 = species[0];
        species1 = species[1];
    }

    /**
     * Configures iterator to return molecules from the set species in the given
     * box.
     */
    public void setBox(Box newBox) {
        if (newBox == null) {
            throw new IllegalArgumentException("You shouldn't pass a null Box.  Why would you do that?");
        }
        box = newBox;
        if (species0.getIndex() > species1.getIndex()) {
            // species were out of order.  swap them
            ISpecies tempSpecies = species0;
            species0 = species1;
            species1 = tempSpecies;
        }
    }

    /**
     * Indicates allowed direction for iteration, relative to specified target
     * atom. If the specified direction is consisent with the direction from the
     * target species to the non-target species (as given by their species index --
     * UP is direction from smaller index to larger index) direction, iteration
     * is performed; if specified direction contradicts species direction, no
     * iteration is performed. Specification of a null direction indicates no
     * limitation, and iteration will be performed if a legitimate target atom
     * is specified.
     */
    public void setDirection(Direction direction) {
        this.direction = direction;
    }

    /**
     * Sets the target molecule with which all pairs are formed. Molecule is
     * determined from the atom specified by the newTargetAtom (which must not
     * be null), which may be the molecule itself or an atom that is part
     * of it. If the atom is not in one of the species given at
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

    /**
     * Finds target molecule as indicated by the target atom. Sets target
     * molecule to null if target atom is null, box is null, or atom is not
     * part of either species.
     */
    private void identifyTargetMolecule() {
        // exception here if you didn't set a target.  call setTarget
        if (targetMolecule.getType() == species0) {
            //target is species0
            allowedDirection = IteratorDirective.Direction.UP;
            iterator = apiUp;
            aiInner.setList(box.getMoleculeList(species1));
        }
        else if (targetMolecule.getType() == species1) {
            //target is species1
            allowedDirection = IteratorDirective.Direction.DOWN;
            iterator = apiDown;
            aiInner.setList(box.getMoleculeList(species0));
        }
        else {
            targetMolecule = null;
        }

        if (direction == null || direction == allowedDirection) {
            aiOuter.setMolecule(targetMolecule);//targetMolecule may be null here
        } else {
            aiOuter.setMolecule(null);
        }
    }

    public int nBody() {
        return 2;
    }
    
    public IMoleculeList next() {
        return iterator.next();
    }
    
    public void reset() {
        identifyTargetMolecule();
        iterator.reset();
    }

    public int size() {
        return iterator.size();
    }
    
    public void unset() {
        iterator.unset();
    }
    
    private static final long serialVersionUID = 1L;
    private final MoleculeIteratorArrayListSimple aiInner;
    private final MoleculeIteratorSinglet aiOuter;
    private ISpecies species0, species1;
    private final MpiInnerFixed apiUp, apiDown;
    private MpiInnerFixed iterator;
    private IteratorDirective.Direction direction, allowedDirection;
    private Box box;
    private IMolecule targetMolecule;
}
