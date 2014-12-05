/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.potential.PotentialMolecular;
import etomica.space.ISpace;

/**
 * @author David Kofke
 *
 * Pair potential given according to the Mayer bonds in a cluster integral.
 * Does not require that the value of the cluster is non-negative.
 */
public class P0Cluster extends PotentialMolecular {

    private static final long serialVersionUID = 1L;
    private BoxCluster boxCluster;
	/**
	 * Constructor for P0Cluster.
	 */
	public P0Cluster(ISpace space) {
		super(0,space);
	}
	
    // let's all pretend that the cluster weight is the energy.
	public double energy(IMoleculeList atoms) {
        return 0;
	}

    public double weight() {
        return boxCluster.getSampleCluster().value(boxCluster);
    }

    public void setBox(IBox box) {
    	boxCluster = (BoxCluster)box;
    }

    public double getRange() {
        return 0;
    }
}
