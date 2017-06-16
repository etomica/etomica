/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;

import etomica.space.Vector;

/**
 * Polytope having parallel faces, and which can be specified via a set of edge vectors.
 */
public interface Parallelotope {

    /**
     * Specifies the polytope via an array of edge vectors.  The number of
     * vectors should be equal to the dimension of the polytope.
     */
    public void setEdgeVectors(Vector[] edgeVectors);
}
