package etomica.virial;

import etomica.atom.iterator.AtomsetIterator;
import etomica.api.IPotential;
import etomica.potential.PotentialCalculation;

/**
 * Calculates the cluster weight associated with current configuration.
 */

public class PotentialCalculationClusterWeightSum extends PotentialCalculation {
        
    public void reset() {
        weight = 1.0;
    }
    public double weight() {
        return weight;
    }

    protected void doCalculation(AtomsetIterator iterator, IPotential potential) {
        if (potential instanceof P0Cluster) {
            weight *= ((P0Cluster)potential).weight();
        }
	}

    protected double weight;
}
