/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;

import etomica.api.IVectorMutable;
import etomica.space.ISpace;

/**
 * Single-atom basis with the coordinate at the origin.
 */
public class BasisMonatomic extends Basis {
    
    /**
     * Creates a single-atom basis with the coordinate at the origin.
     */
    public BasisMonatomic(ISpace space) {
        super(new IVectorMutable[] {space.makeVector()});
    }
    	
    private static final long serialVersionUID = 1L;
}
