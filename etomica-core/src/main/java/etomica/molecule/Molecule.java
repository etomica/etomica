/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.species.ISpecies;

public class Molecule implements IMolecule, java.io.Serializable {

    public Molecule(ISpecies species, int numLeafAtoms) {
        this.species = species;
        childList = new AtomArrayList(numLeafAtoms);
    }
    
    /**
     * Returns a string of digits that uniquely identifies this atom.  String is
     * formed by concatenating the ordinal of this atom to the signature
     * given by the parent of this atom.  If atom has no parent, forms a string
     * from only the ordinal.
     */
    public String signature() {
        return species.getIndex()+" "+index;
    }

    /**
     * Returns a string formed by concatenating the signature of this atom
     * to a string that identifies it as a molecule.
     */
    public final String toString() {
        return "Molecule(" + getType().getIndex()+" "+index + ")";
    }

    /**
     * Adds the given Atom as a child of this Atom.  The given child Atom
     * should be parentless when this method is called.
     * @throws IllegalArgumentException if the given atom already has a parent.
     */
    public void addChildAtom(IAtom newChildAtom) {
        if(newChildAtom.getParentGroup() != null) {//new parent is null
            throw new IllegalArgumentException(newChildAtom+" is already the child of "+newChildAtom.getParentGroup());
        }

        newChildAtom.setParent(this);

        newChildAtom.setIndex(childList.size());
        childList.add(newChildAtom);
    }
    
    /**
     * Removes the given child Atom from this AtomGroup.
     * @throws IllegalArgumentException if the given atom is not a child.
     */
    public void removeChildAtom(IAtom oldChildAtom) {
        for (int i = 0; i<childList.size(); i++) {
            if (childList.get(i) == oldChildAtom) {
                oldChildAtom.setParent(null);
                childList.removeAndReplace(i);
                childList.maybeTrimToSize();
                if (childList.size() > i) {
                    // reassign the old last Atom (which is now in the removed
                    // Atom's place) to have the old Atom's index.
                    childList.get(i).setIndex(i);
                }
                return;
            }
        }
        throw new IllegalArgumentException(oldChildAtom+" is not a child");
    }

    /**
     * @return the childList
     */
    public final IAtomList getChildList() {
        return childList;
    }
    
    public final void setIndex(int newIndex) {
        index = newIndex;
    }
    
    public final int getIndex() {
        return index;
    }
    
    public final ISpecies getType() {
        return species;
    }

    private static final long serialVersionUID = 1L;
    
    protected int index;
    protected final AtomArrayList childList;
    protected final ISpecies species;
}
