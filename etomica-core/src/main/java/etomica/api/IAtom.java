/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

/**
 * IAtom is an interface for atoms in the simulation.  Atoms have a type
 * (IAtomType), a position (IVectorMutable) and a parent molecule (IMolecule).
 * Atoms also have indices (one for the list of atoms in the box, one for the
 * list of atoms in a molecule).
 */
public interface IAtom {

    /**
     * @return this IAtom's index, which is its place in the parent AtomGroup's
     * list of child IAtoms.
     */
    public int getIndex();

    /**
     * Informs the IAtom of its index within the atoms of its parent molecule.
     * This should only be called by the parent molecule.
     * @param index the new index
     */
    public void setIndex(int index);

    /**
     * Sets the atom's global index to the give value.  This method should only
     * be called by the IBox.
     */
    public void setLeafIndex(int newGlobalIndex);

    /**
     * Returns the global index (within the Box) of this Atom.  The global
     * index is unique to the IAtom in the Box.  The IAtom's global index may
     * change over the course of a simulation due to addition or removal of
     * other IAtoms in the Box.  An BoxGlobalAtomIndexEvent is fired by
     * the Box's event manager when an Atom's global index changes. 
     */
    public int getLeafIndex();

    /**
     * Informs the Atom that the given IMolecule is its parent.
     * This method should only be called by the parent molecule.
     * @param newParent the new parent molecule
     */
    public void setParent(IMolecule newParent);

    /**
     * @return the parent molecule of this IAtom.
     */
    public IMolecule getParentGroup();

    /**
     * @return the Atom type, holding properties held in common with other 
     * atoms of the same type.
     */
    public IAtomType getType();

    /**
     * @return the position of the IAtom.  Modifying the returned vector will
     * alter the IAtom's position.
     */
    public IVectorMutable getPosition();
}
