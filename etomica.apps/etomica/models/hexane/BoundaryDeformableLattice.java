/*
 * Created on Jun 15, 2005
 */
package etomica.models.hexane;

import etomica.lattice.crystal.Primitive;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Vector;

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
        super(primitive.space, makeDimensionVectors(primitive, nCells));
    }
    
    /**
     * Creates boundary vectors in the direction of the primitive vectors, 
     * scaled by nCells
     */
    private static final Vector[] makeDimensionVectors(Primitive primitive, int[] numCells) {
        Vector[] dimensionVectors = new Vector[primitive.space.D()];
        for (int i=0; i<dimensionVectors.length; i++) {
            dimensionVectors[i] = primitive.space.makeVector();
            dimensionVectors[i].Ea1Tv1(numCells[i],primitive.vectors()[i]);
        }
        return dimensionVectors;
    }

    private static final long serialVersionUID = 1L;
}
