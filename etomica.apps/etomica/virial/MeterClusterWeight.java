package etomica.virial;

import etomica.DataInfo;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.data.DataSourceScalar;
import etomica.units.Dimension;

/**
 * Meter to calculate the sampling weight of a cluster configuration.   
 */
 
public class MeterClusterWeight extends DataSourceScalar {
    
    public MeterClusterWeight(PotentialMaster potentialMaster) {
        super(new DataInfo("Cluster Weight",Dimension.NULL));
        potential = potentialMaster;
    }
      
    public Dimension getDimension() {return Dimension.NULL;}
    
    public double getDataAsScalar() {
    	weight.reset();
    	potential.calculate(phase, up, weight);
    	return weight.sum();
    }

    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    private Phase phase;
    private final PotentialMaster potential;
    private final PotentialCalculationClusterWeightSum weight = new PotentialCalculationClusterWeightSum();
    private final IteratorDirective up = new IteratorDirective();
}