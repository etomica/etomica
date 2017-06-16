/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.space;

import etomica.lattice.crystal.Primitive;

/**
 * Deformable boundary that takes the shape of a primitive with some number
 * of cells in each dimension.
 */
public class BoundaryDeformableLattice extends BoundaryDeformablePeriodic {
    
    /**
     * Creates boundary with the shape of the given primitive and the given
     * number of cells for each dimension.  Changes to the primitive after
     * creating this instance will not affect the boundary.
     */
    public BoundaryDeformableLattice(Primitive primitive, int[] nCells){
        super(primitive.getSpace(), makeDimensionVectors(primitive, nCells));
    }

    /**
     * Creates boundary with the shape of the given primitive and the given
     * number of cells for each dimension.  Changes to the primitive after
     * creating this instance will not affect the boundary.
     */
    public BoundaryDeformableLattice(Primitive primitive, double[] nCells){
        super(primitive.getSpace(), makeDimensionVectors(primitive, nCells));
    }

    /**
     * Creates boundary vectors in the direction of the primitive vectors, 
     * scaled by nCells
     */
    private static final Vector[] makeDimensionVectors(Primitive primitive, int[] numCells) {
        Vector[] dimensionVectors = new Vector[primitive.getSpace().D()];
        for (int i=0; i<dimensionVectors.length; i++) {
            dimensionVectors[i] = primitive.getSpace().makeVector();
            dimensionVectors[i].Ea1Tv1(numCells[i],primitive.vectors()[i]);
        }
        return dimensionVectors;
    }

    /**
     * Creates boundary vectors in the direction of the primitive vectors, 
     * scaled by nCells
     */
    private static final Vector[] makeDimensionVectors(Primitive primitive, double[] numCells) {
        Vector[] dimensionVectors = new Vector[primitive.getSpace().D()];
        for (int i=0; i<dimensionVectors.length; i++) {
            dimensionVectors[i] = primitive.getSpace().makeVector();
            dimensionVectors[i].Ea1Tv1(numCells[i],primitive.vectors()[i]);
//          dimensionVectors[i].E(numCells[i]);
        }
        return dimensionVectors;
    }

    private static final long serialVersionUID = 1L;
}
