/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.threaded.domain;

import etomica.box.Box;
import etomica.nbr.list.NeighborListManager;
import etomica.space.Space;

public class NeighborListManagerThreaded extends NeighborListManager {

    public NeighborListManagerThreaded(PotentialMasterListThreaded potentialMasterListThreaded,
                                       double range, Box box, Space space) {
        super(potentialMasterListThreaded, range, box, space);
        
    }

    
    protected void neighborSetup() {
      
        NeighborCellManagerThreaded cellManager = ((PotentialMasterListThreaded)potentialMaster).getNbrCellManagerThreaded(box);
        cellManager.assignCellAll();
        super.neighborSetup();
    }
   



}
