/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.api.ISpecies;
import etomica.atom.MoleculeArrayList;
import etomica.atom.iterator.IteratorDirective.Direction;

/**
 * Iterator for the molecules of a single species in a box.  Can be targeted to
 * give the molecule containing a specific atom (if consistent with the species),
 * or all molecules of the species in a box.<br>
 * This class is used by PotentialMaster to iterate over molecules for single-body
 * potentials.
 */
public class MoleculeIteratorMolecule extends MoleculeIteratorArrayListSimple implements
        MoleculesetIteratorPDT, MoleculeIteratorBoxDependent {

    /**
     * @param species species for which molecules are returned as iterates.
     * species must not be null.
     */
    public MoleculeIteratorMolecule(ISpecies species) {
        super();
        this.species = species;
    }

    /**
     * Sets the box containing the molecules for iteration. A null
     * box conditions iterator to give no iterates.
     * @throws NullPointerException if the Box is null
     */
    public void setBox(Box newBox) {
        box = newBox;
    }
    
    public void reset() {
        setList();
        super.reset();
    }

    /**
     * Specifies molecule returned by this iterator, as the one containing
     * the given target atom.  Only the first element of the array is relevant.
     * If argument is null, of zero length, or if targetAtom[0] is null, then
     * no target atom is specified and all molecules of the species in the
     * box will be given on iteration.
     * 
     * @throws NullPointerException
     *          if targetAtom is null
     * @throws IllegalArgumentException
     *          if targetAtom.count() is not 0 or 1
     */
    public void setTarget(IMolecule newTargetAtom) {
        targetAtom = newTargetAtom;
    }

    /** 
     * Has no effect, but is included as part of the AtomsetIteratorPDT interface.
     */
    public void setDirection(Direction direction) {
        //ignore
    }

    /**
     * Configures the list iterator with a list appropriate to the specified
     * box and target.
     */
    private void setList() {
        if(targetAtom == null) {
            setList(box.getMoleculeList(species));
        
        //target specified -- give it as only iterate if descended from species
        } else {
            littleList.clear();
            if (targetAtom.getType() == species) {
                littleList.add(targetAtom);
            }
            setList(littleList);
        }
    }

    private static final long serialVersionUID = 1L;
    private final ISpecies species;
    private final MoleculeArrayList littleList = new MoleculeArrayList();
    private Box box;
    private IMolecule targetAtom;
}
