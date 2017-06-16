/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.potential.PotentialMaster;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.PotentialMasterList;

/**
 * A PotentialGroup class to work with PotentialMasterHybrid.  Holds a reference
 * to cell, list and hybrid PotentialMasters so that all of them can be notified
 *
 * @author Andrew Schultz
 */
public class PotentialGroupHybrid extends PotentialGroupNbr {

    public PotentialGroupHybrid(int nBody) {
        super(nBody);
    }

    public void setPotentialMaster(PotentialMaster newPotentialMaster) {
        // we'll get notification of each PotentialMaster instance (cell, list, 
        // hybrid), but we only care about the hyrbid one.  It will notify the
        // others.
        if (newPotentialMaster instanceof PotentialMasterHybrid) {
            // we got what we wanted.  
            super.setPotentialMaster(newPotentialMaster);
        }
    }
    
    private static final long serialVersionUID = 1L;
    protected PotentialMasterCell potentialMasterCell;
    protected PotentialMasterList potentialMasterList;
}
