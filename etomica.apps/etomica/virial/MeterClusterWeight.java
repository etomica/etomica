package etomica.virial;

import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataSourceScalar;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.units.Dimension;
import etomica.units.Null;

/**
 * Meter to calculate the sampling weight of a cluster configuration.   
 */
 
public class MeterClusterWeight extends DataSourceScalar {
    
    public MeterClusterWeight(PotentialMaster potentialMaster) {
        super("Cluster Weight",Null.DIMENSION);
        potential = potentialMaster;
    }
      
    public Dimension getDimension() {return Null.DIMENSION;}
    
    public double getDataAsScalar() {
    	weight.reset();
    	potential.calculate(phase, up, weight);
    	return weight.weight();
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