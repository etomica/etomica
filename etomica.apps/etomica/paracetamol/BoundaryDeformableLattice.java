
package etomica.paracetamol;

import etomica.lattice.crystal.Primitive;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.IVector;
import etomica.util.IRandom;

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
    public BoundaryDeformableLattice(Primitive primitive, IRandom random, int[] nCells){
        super(primitive.getSpace(), random, makeDimensionVectors(primitive, nCells));
    }
    
    /**
     * Creates boundary vectors in the direction of the primitive vectors, 
     * scaled by nCells
     */
    private static final IVector[] makeDimensionVectors(Primitive primitive, int[] numCells) {
        IVector[] dimensionVectors = new IVector[primitive.getSpace().D()];
        for (int i=0; i<dimensionVectors.length; i++) {
            dimensionVectors[i] = primitive.getSpace().makeVector();
            dimensionVectors[i].Ea1Tv1(numCells[i],primitive.vectors()[i]);
        }
        return dimensionVectors;
    }

    private static final long serialVersionUID = 1L;
}
