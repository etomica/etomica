/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
