/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential; 

import etomica.nbr.NeighborCriterion;
import etomica.space.Space;

/**
 * Potential that does not depend on any atom positions.
 * Typically used to implement long-range corrections for potential truncation.
 * Potential thus depends on box parameters, such as the number of molecules 
 * and the volume.
 *
 * @author David Kofke
 */

public abstract class Potential0 extends Potential {
      
    public Potential0(Space space) {
        super(0, space);
    }
    
    /**
     * Returns zero.
     */
    public double getRange() {
        return 0.0;
    }
    
    public NeighborCriterion getCriterion() {
        return null;
    }
                        
}//end of Potential0



