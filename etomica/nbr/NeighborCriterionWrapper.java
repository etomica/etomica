package etomica.nbr;

import etomica.Phase;


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
    
    private final NeighborCriterion[] criteria;
    private double neighborRange;
}
