/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.potential.PotentialMaster;
import etomica.data.meter.MeterPotentialEnergy;

public class ClusterSumExternalField extends ClusterSum {

	public ClusterSumExternalField(ClusterBonds[] subClusters,
			double[] subClusterWeights, MayerFunction[] fArray) {
		super(subClusters, subClusterWeights, fArray);
		// TODO Auto-generated constructor stub
	}
	public void setPotentialMaster(PotentialMaster potentialMaster){
		meterPE=new MeterPotentialEnergy(potentialMaster);
	}
	public double value(BoxCluster box) {
		meterPE.setBox(box);
		return super.value(box);
	}
	protected void calcValue() {
		super.calcValue();
	
		double u=meterPE.getDataAsScalar();
		value*=Math.exp(-u*beta)-1;
	}
	protected MeterPotentialEnergy meterPE; 
}
