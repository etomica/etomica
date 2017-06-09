/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.box.Box;
import etomica.atom.IMoleculeList;
import etomica.api.ISpecies;
import etomica.atom.MoleculePair;

/**
 * Gives pairs formed from the molecules of two different species in a box.
 * Species are specified at construction and cannot be changed afterwards.
 */

public class MpiInterspeciesAA implements MoleculesetIteratorBoxDependent {

    /**
     * Constructs iterator that provides iterates taken from the molecules of
     * two species. Given array is sorted in increasing order of species index. 
     * Then atoms in iterate pairs will be such that species of atom0 is
     * species[0] (having smaller species index), and species of atom1 is species[1]
     * (having larger species index).
     * 
     * @param species
     *            array of two different, non-null species
     * @throws NullPointerException
     *             if species or one of its elements is null
     * @throws IllegalArgumentException
     *             if species.length != 2 or if species[0] == species[1]
     */
    public MpiInterspeciesAA(ISpecies[] species) {
        if(species.length != 2) {
            throw new IllegalArgumentException("Incorrect array length; must be 2 but length is "+species.length);
        }

        // we need to sort these.  we'll do that once we have the box
        species0 = species[0];
        species1 = species[1];
        if (species0 == null || species1 == null) {
            throw new NullPointerException(
                    "Constructor of ApiInterspeciesAA requires two non-null species");
        }
        if (species0 == species1) {
            throw new IllegalArgumentException(
                    "Constructor of ApiInterspeciesAA requires two different species");
        }
        
        unset();
    }

    /**
     * Configures iterator to return molecules from the set species in the given
     * box.
     * @throws NullPointerException if the Box is null
     */
    public void setBox(Box box) {
        if (species0.getIndex() > species1.getIndex()) {
            // species were out of order.  swap them
            ISpecies tempSpecies = species0;
            species0 = species1;
            species1 = tempSpecies;
        }
        outerList = box.getMoleculeList(species0);
        innerList = box.getMoleculeList(species1);
        unset();
    }

    /**
     * Sets iterator in condition to begin iteration.
     * 
     * @throws IllegalStateException
     *             if outer and inner lists have been set to the same instance
     */
    public void reset() {
        if (outerList == innerList) {
            throw new IllegalStateException(
                    "ApiInterList will not work correctly if inner and outer lists are the same instance");
        }
        outerIndex = 0;
        if (outerList.getMoleculeCount() == 0) {
            innerIndex = innerList.getMoleculeCount() - 1;
            return;
        }
        innerIndex = -1;
        atoms.atom0 = outerList.getMolecule(outerIndex);
    }

    /**
     * Sets iterator such that hasNext is false.
     */
    public void unset() {
        if (outerList != null) {
            outerIndex = outerList.getMoleculeCount() - 1;
        }
        if (innerList != null) {
            innerIndex = innerList.getMoleculeCount() - 1;
        }
    }

    /**
     * Returns the next iterate pair. Returns null if there are no more
     * iterates.
     */
    public IMoleculeList next() {
        if (innerIndex > innerList.getMoleculeCount() - 2) {
            if (outerIndex > outerList.getMoleculeCount() - 2 || innerList.getMoleculeCount() == 0) {
                return null;
            }
            outerIndex++;
            atoms.atom0 = outerList.getMolecule(outerIndex);
            innerIndex = -1;
        }
        innerIndex++;
        atoms.atom1 = innerList.getMolecule(innerIndex);
        return atoms;
    }

    /**
     * Returns the number of iterates, which is list.size*(list.size-1)/2
     */
    public int size() {
        return outerList.getMoleculeCount() * innerList.getMoleculeCount();
    }

    /**
     * Returns 2, indicating that this is a pair iterator
     */
    public int nBody() {
        return 2;
    }

    private static final long serialVersionUID = 1L;
    private IMoleculeList outerList, innerList;
    private int outerIndex, innerIndex;
    private final MoleculePair atoms = new MoleculePair();
    private ISpecies species0, species1;
}
