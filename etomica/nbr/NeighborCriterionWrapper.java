package etomica.nbr;

import etomica.NearestImageVectorSource;
import etomica.Phase;
import etomica.nbr.cell.AtomsetIteratorCellular;


/**
 * Fake criterion to just hold a range to give to PotentialMasterNbr
 */

public class NeighborCriterionWrapper extends NeighborCriterionAll {

    public NeighborCriterionWrapper(NeighborCriterion[] subCriteria) {
        criteria = subCriteria;
    }
    
    public void setNeighborRange(double range) {
        neighborRange = range;
    }
    
    public double getNeighborRange() {
        return neighborRange;
    }
    
    public void setPhase(Phase p) {
        for (int i=0; i<criteria.length; i++) {
            criteria[i].setPhase(p);
        }
    }
    
    public void setNearestImageVectorSource(NearestImageVectorSource nivs) {
        for (int i=0; i<criteria.length; i++) {
            criteria[i].setNearestImageVectorSource(nivs);
        }
    }
    
    public void setCellIterator(AtomsetIteratorCellular api) {
        for (int i=0; i<criteria.length; i++) {
            criteria[i].setCellIterator(api);
        }
    }

    private final NeighborCriterion[] criteria;
    private double neighborRange;
}
