package etomica.math.geometry;

import etomica.api.IVector;

/**
 * Polytope having parallel faces, and which can be specified via a set of edge vectors.
 */
public interface Parallelotope {

    /**
     * Specifies the polytope via an array of edge vectors.  The number of
     * vectors should be equal to the dimension of the polytope.
     */
    public void setEdgeVectors(IVector[] edgeVectors);
}
