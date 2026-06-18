/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.Debug;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;

 /**
  * Atom that represents a physical atom with a position.
  * <p>
  * @author David Kofke, Andrew Schultz, and C. Daniel Barnes
  * 
  */
public class Atom implements IAtom {

     protected final AtomType type;
     protected final Vector position;
     protected int index;
     protected IMolecule parent;
     protected int leafIndex;

     public Atom(Space space, AtomType type) {
        super();
        this.type = type;
        position = space.makeVector();
    }

    /**
     * Makes a simple atom for the given space.  Coordinate is non-kinetic sphere;
     * node is for a leaf atom; type is a sphere with unit mass and unit size,
     * unique to the new atom; depth is 0.
     */
    public Atom(Space space) {
        super();
        type = null;
        position = space.makeVector();
    }

    /**
     * Returns a string of digits that uniquely identifies this atom.  String is
     * formed by concatenating the ordinal of this atom to the signature
     * given by the parent of this atom.  If atom has no parent, forms a string
     * from only the ordinal.
     */
    public String signature() {
        if(parent != null) {
            if(parent instanceof Molecule) {
                return ((Molecule)parent).signature() + " " + getIndex();
            }
        }
        return Integer.toString(getIndex());
    }

    /**
     * Returns a string formed by concatenating the signature of this atom
     * to a string that identifies it as a species master, species agent,
     * molecule, group, or (leaf) atom.
     */
    public final String toString() {
        return "Atom(" + signature() + ")";
    }

    public final int getIndex() {
        return index;
    }

    public final void setIndex(int newIndex) {
        index = newIndex;
    }

    /**
     * Informs the Atom that the given AtomGroup is its parent.
     * This method should only be called by the parent.
     */
    public void setParent(IMolecule newParent) {
        if (Debug.ON && ((newParent == null) == (parent == null))) {
            // newParent must be null (removal) if current parent is not null
            // new Parent must not be null if current parent null
            throw new IllegalArgumentException(newParent+" is not my parent");
        }
        parent = newParent;
    }

    public IMolecule getParentGroup() {
        return parent;
    }

    public Vector getPosition() {
        return position;
    }

    public int getLeafIndex() {
        return leafIndex;
    }

    public void setLeafIndex(int newLeafIndex) {
         leafIndex = newLeafIndex;
     }

    /**
     * @return the Atom type, holding properties held in common with other
     * atoms made by this atom's factory.
     */
    public final AtomType getType() {
        return type;
    }

    public void copyCoordinatesFrom(IAtom atom) {
        position.E(atom.getPosition());
    }

    @Override
    public void saveState(Writer fw) throws IOException {
        int D = position.getD();
        fw.write("" + leafIndex);
        for (int i = 0; i < D; i++) fw.write(" " + position.getX(i));
        fw.write("\n");
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        String[] bits = br.readLine().split(" ");
        leafIndex = Integer.parseInt(bits[0]);
        int D = position.getD();
        for (int i = 0; i < D; i++) position.setX(i, Double.parseDouble(bits[1 + i]));
    }
}
