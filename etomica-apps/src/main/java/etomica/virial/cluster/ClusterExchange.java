/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import etomica.potential.IPotentialField;
import etomica.virial.BoxCluster;


public class ClusterExchange implements ClusterAbstract {
	protected double beta;
	protected IPotentialField p1;
	public ClusterExchange (IPotentialField p1) {
		this.p1 = p1;		
	}

	public ClusterAbstract makeCopy() {
		ClusterExchange cEx = new ClusterExchange(p1);
		cEx.setTemperature(1/beta);
		return cEx;
	}

	public int pointCount() {
		return 1;
	}

	public double value(BoxCluster box) {
		int beads = box.getLeafList().size();
		double sum = 0;
		for (int i=0; i<beads; i++) {
			sum += p1.u(box.getLeafList().get(i));
		}
		sum /= beads;
		sum = Math.exp(-beta*sum);
		return sum;
	}

	public void setTemperature(double temperature) {
		beta = 1/temperature;
	}

}
