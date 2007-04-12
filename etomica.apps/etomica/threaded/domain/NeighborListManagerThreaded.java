package etomica.threaded.domain;

import etomica.nbr.list.NeighborListManager;
import etomica.phase.Phase;

public class NeighborListManagerThreaded extends NeighborListManager {

    public NeighborListManagerThreaded(PotentialMasterListThreaded potentialMasterListThreaded,
            double range, Phase phase) {
        super(potentialMasterListThreaded, range, phase);
        
    }

    
    protected void neighborSetup() {
      
        NeighborCellManagerThreaded cellManager = ((PotentialMasterListThreaded)potentialMaster).getNbrCellManagerThreaded(phase);
        cellManager.assignCellAll();
        super.neighborSetup();
    }
   



}
