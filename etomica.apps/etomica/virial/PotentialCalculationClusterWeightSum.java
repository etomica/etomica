package etomica.virial;

import etomica.AtomsetIterator;
import etomica.Potential;
import etomica.potential.PotentialCalculation;

/**
 * Calculates the cluster weight associated with current configuration.
 */

public class PotentialCalculationClusterWeightSum extends PotentialCalculation {
    protected double sum = 0.0;
        
    public void reset() {sum = 1.0;}
    public double sum() {return sum;}

    protected void doCalculation(AtomsetIterator iterator, Potential potential) {
        sum = ((P0Cluster)potential).weight();
	}

}
