package etomica.virial;

import etomica.Phase;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMaster;
import etomica.units.Dimension;

/**
 * Meter to calculate the sampling weight of a cluster configuration.   
 */
 
public class MeterClusterWeight extends DataSourceScalar {
    
    public MeterClusterWeight(PotentialMaster potentialMaster) {
        super("Cluster Weight",Dimension.NULL);
        potential = potentialMaster;
    }
      
    public Dimension getDimension() {return Dimension.NULL;}
    
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