package etomica.threaded.domain;

import etomica.nbr.list.NeighborListManager;
import etomica.box.Box;

public class NeighborListManagerThreaded extends NeighborListManager {

    public NeighborListManagerThreaded(PotentialMasterListThreaded potentialMasterListThreaded,
            double range, Box box) {
        super(potentialMasterListThreaded, range, box);
        
    }

    
    protected void neighborSetup() {
      
        NeighborCellManagerThreaded cellManager = ((PotentialMasterListThreaded)potentialMaster).getNbrCellManagerThreaded(box);
        cellManager.assignCellAll();
        super.neighborSetup();
    }
   



}
