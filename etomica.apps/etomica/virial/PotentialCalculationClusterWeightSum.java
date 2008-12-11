package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;
import etomica.potential.PotentialCalculationMolecular;

/**
 * Calculates the cluster weight associated with current configuration.
 */

public class PotentialCalculationClusterWeightSum implements PotentialCalculationMolecular {
        
    public void reset() {
        weight = 1.0;
    }
    public double weight() {
        return weight;
    }

    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
    }
    
    public void doCalculation(IMoleculeList atoms, IPotentialMolecular potential) {
        // we'll also get intramolecular potentials... ignore them
        if (potential instanceof P0Cluster) {
            weight *= ((P0Cluster)potential).weight();
        }
	}

    protected double weight;
}
