package etomica.nbr;


/**
 * Fake criterion to just hold a range to give to PotentialMasterNbr
 */

public class NeighborCriterionFake extends NeighborCriterionAll {

    public void setNeighborRange(double range) {
        neighborRange = range;
    }
    
    public double getNeighborRange() {
        return neighborRange;
    }
    
    private double neighborRange;
}
