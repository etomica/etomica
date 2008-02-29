package etomica.threaded.domain;

import etomica.nbr.list.NeighborListManager;
import etomica.box.Box;
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
