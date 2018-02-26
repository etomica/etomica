/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

import etomica.space.Vector;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;

/**
 * A lattice with sites given by the "atom" sites of a crystal. Sites of this
 * lattice are instances of Vector. The dimension of a BravaisLatticeCrystal is
 * one more than the dimension of the underlying Bravais lattice forming the
 * crystal; the extra index specifies the basis atom at the site referenced by
 * the other indices.
 */
public class BravaisLatticeCrystal extends BravaisLattice {

    /**
     * Constructs a lattice having sites given by the "atom" sites of the given
     * crystal. Note: the positions defined by the Basis are specified in the
     * coordinate system of the given primitive, which typically will not be the
     * rectangular cartesian coordinates. This will cause the basis
     * to in effect be rescaled when the primitive is resized.
     */
    public BravaisLatticeCrystal(Primitive primitive, Basis basis) {
        super(primitive);
        this.basis = basis;
        D = primitive.getSpace().D() + 1;
        crystalIndex = new int[D - 1];
        position = getSpace().makeVector();
    }

    /**
     * Returns the spatial dimension + 1. The extra index specifies the basis
     * atom.
     */
    public int D() {
        return D;
    }

    /**
     * Returns a Vector instance giving the location of the referenced site. The
     * first D-1 indices indicate the Bravais-lattice position, and the last
     * index specifies the basis atom at the Bravais site. The same Vector
     * instance is returned with each call.
     */
    public Vector site(int[] index) {
        if (index.length != D)
            throw new IllegalArgumentException(
                    "index given to site method of lattice must have number of elements equal to dimension of lattice");
        System.arraycopy(index, 0, crystalIndex, 0, D - 1);
        Vector latticePosition = (Vector) super.site(crystalIndex);
        Vector basisCoordinate = basis.getScaledCoordinates()[index[D - 1]];
        Vector[] primitiveVectors = primitive.vectors();
        position.E(latticePosition);
        for (int i = 0; i < basisCoordinate.getD(); i++) {
            position.PEa1Tv1(basisCoordinate.getX(i), primitiveVectors[i]);// basis is specified in the frame defined by the primitive vectors
        }
        return position;
    }

    public Basis getBasis() {
        return basis;
    }

    private static final long serialVersionUID = 1L;
    protected final Basis basis;
    private final int[] crystalIndex;
    private final int D;
    private final Vector position;
}
