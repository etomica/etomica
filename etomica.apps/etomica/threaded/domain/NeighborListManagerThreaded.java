package etomica.threaded.domain;

import etomica.nbr.list.NeighborListManager;
import etomica.api.IBox;
import etomica.space.ISpace;

public class NeighborListManagerThreaded extends NeighborListManager {

    public NeighborListManagerThreaded(PotentialMasterListThreaded potentialMasterListThreaded,
            double range, IBox box, ISpace space) {
        super(potentialMasterListThreaded, range, box, space);
        
    }

    
    protected void neighborSetup() {
      
        NeighborCellManagerThreaded cellManager = ((PotentialMasterListThreaded)potentialMaster).getNbrCellManagerThreaded(box);
        cellManager.assignCellAll();
        super.neighborSetup();
    }
   



}
