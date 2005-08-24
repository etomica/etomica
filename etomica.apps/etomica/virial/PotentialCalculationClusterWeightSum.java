package etomica.virial;

import etomica.Potential;
import etomica.atom.AtomsetIterator;
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

    protected void doCalculation(AtomsetIterator iterator, Potential potential) {
        iterator.reset();
        while (iterator.hasNext()) {
            weight *= potential.energy(iterator.next());
        }
	}

    protected double weight;
}
