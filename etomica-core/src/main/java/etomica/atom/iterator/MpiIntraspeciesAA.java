/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import java.io.Serializable;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.atom.MoleculePair;
import etomica.util.Debug;

/**
 * Gives pairs formed from the molecules of a single species in a box. Species
 * is specified at construction and cannot be changed afterwards.
 */
public class MpiIntraspeciesAA implements MoleculesetIteratorBoxDependent, Serializable {

    /**
     * @param species
     *            species whose molecules will form the pair iterates
     */
    public MpiIntraspeciesAA(ISpecies species) {
        if (species == null) {
            throw new NullPointerException("Constructor of ApiIntraspecies1A requires a non-null species");
        }
        this.species = species;
    }

    /**
     * Configures iterator to return molecules from the set species in the given
     * box.
     * @throws NullPointerException if the Box is null
     */
    public void setBox(IBox box) {
        list = box.getMoleculeList(species);
        unset();
    }

    /**
     * Sets iterator in condition to begin iteration.
     */
    public void reset() {
        if (list.getMoleculeCount() < 2) {
            outerIndex = 2;
            innerIndex = 2;
            return;
        }
        outerIndex = 0;
        innerIndex = 0;
        atoms.atom0 = list.getMolecule(0);
    }

    /**
     * Sets iterator such that next is null.
     */
    public void unset() {
        outerIndex = list.getMoleculeCount() - 2;
        innerIndex = list.getMoleculeCount() - 1;
    }

    /**
     * Returns the number of iterates, which is list.size*(list.size-1)/2
     */
    public int size() {
        return list.getMoleculeCount() * (list.getMoleculeCount() - 1) / 2;
    }

    /**
     * Returns the next iterate pair. Returns null if hasNext() is false.
     */
    public IMoleculeList next() {
        if (innerIndex > list.getMoleculeCount() - 2) {
            if (outerIndex > list.getMoleculeCount() - 3) {
                return null;
            }
            outerIndex++;
            atoms.atom0 = list.getMolecule(outerIndex);
            innerIndex = outerIndex;
        }
        innerIndex++;
        atoms.atom1 = list.getMolecule(innerIndex);
        if (Debug.ON && atoms.atom0 == atoms.atom1) {
            throw new RuntimeException("oops");
        }
        return atoms;
    }

    /**
     * Returns 2, indicating that this is a pair iterator
     */
    public int nBody() {
        return 2;
    }

    private static final long serialVersionUID = 1L;
    private IMoleculeList list;
    private int outerIndex, innerIndex;
    private final MoleculePair atoms = new MoleculePair();
    private final ISpecies species;
}
