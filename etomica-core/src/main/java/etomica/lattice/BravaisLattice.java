/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

import etomica.space.Vector;
import etomica.lattice.crystal.Primitive;
import etomica.space.Space;

/**
 * Arbitrary-dimension Bravais Lattice, in which the sites are instances of 
 * etomica.space.Vector, with positions given as linear combinations of a set of
 * primitive vectors.
 */

public class BravaisLattice implements SpaceLattice, java.io.Serializable {

    public BravaisLattice(Primitive primitive) {
        this.primitive = primitive;
        latticeVector = primitive.getSpace().makeVector();
    }

    public int D() {
        return getSpace().D();
    }
    
    public Space getSpace() {
        return primitive.getSpace();
    }
    
    /**
     * Calculates and returns a vector that is the spatial position given
     * by adding together the primitive vectors, each multiplied by the corresponding
     * integer index given by the array argument.  The returned object is an instance of
     * Vector, and the same instance is returned with every call. The returned vector is configured
     * on-the-fly from the the lattice primitive and the given index.  Index may comprise 
     * any integer values (positive, negative, or zero).
     */
    public Object site(int[] index) {
        if(index.length != getSpace().D()) throw new IllegalArgumentException("index given to site method of lattice must have number of elements equal to dimension of lattice");
        latticeVector.E(0);
        Vector[] latticeVectors = primitive.vectors();
        for(int i=0; i<index.length; i++) {
            latticeVector.PEa1Tv1(index[i], latticeVectors[i]);
        }
        return latticeVector;
    }

    /**
     * Sets the primitive for this lattice to the one given, and
     * updates the site positions.
     */
    public void setPrimitive(Primitive primitive) {
        this.primitive = primitive;
    }
    
    /**
     * Returns the primitive object used to construct this lattice. The returned
     * object is not a clone, and changes to it will affect the lattice.
     */
    public Primitive getPrimitive() {
        return primitive;
    }
    
    public double[] getLatticeConstants() {
        return primitive.getSize();
    }

    private static final long serialVersionUID = 1L;
    protected Primitive primitive;
    private final Vector latticeVector;
    
}
