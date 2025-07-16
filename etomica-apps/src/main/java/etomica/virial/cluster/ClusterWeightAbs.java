/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;


import etomica.molecule.CenterOfMass;
import etomica.virial.BoxCluster;

/**
 * @author andrew
 *
 * cluster weight wrapper (absolute value of wrapped cluster)
 */
public class ClusterWeightAbs implements ClusterWeight, java.io.Serializable {
	
    protected final ClusterAbstract weightCluster;
    protected boolean doAbs = true;
	protected double minValue = 0;
	protected double maxDistance = -1;
	
	public ClusterWeightAbs(ClusterAbstract cluster) {
		weightCluster = cluster;
	}

    /**
     * Sets the cluster to actually return the |value| (as default) or just value.
     */
    public void setDoAbs(boolean doAbs) {
        this.doAbs = doAbs;
    }

	/**
	 * Sets a minimum return value (for use with abs).  This can help find an
	 * initial configuration at the start of a sim.
	 */
	public void setMinValue(double minValue) {
		this.minValue = minValue;
	}

	public static ClusterWeight makeWeightCluster(ClusterAbstract cluster) {
		if (cluster instanceof ClusterWeight) {
			return (ClusterWeight)cluster;
		}
		return new ClusterWeightAbs(cluster);
	}

	public static ClusterWeight[] makeWeightClusters(ClusterAbstract[] clusters) {
		ClusterWeight[] weightClusters = new ClusterWeight[clusters.length];
		for (int i=0; i<clusters.length; i++) {
			weightClusters[i] = makeWeightCluster(clusters[i]);
		}
		return weightClusters;
	}
    
    public ClusterAbstract getSubCluster() {
        return weightCluster;
    }
    
    public ClusterAbstract makeCopy() {
        return new ClusterWeightAbs(weightCluster.makeCopy());
    }
	
	public int pointCount() {
		return weightCluster.pointCount();
	}
	
	public double value(BoxCluster box) {
		if (maxDistance >= 0) {
			for (int i=0; i<box.getMoleculeList().size(); i++) {
				double dist = CenterOfMass.position(box, box.getMoleculeList().get(i)).squared();
				if (dist > maxDistance*maxDistance) {
					return 0;
				}
			}

		}
		double v = weightCluster.value(box);
		System.out.println(v);
		return doAbs ? Math.max(minValue, Math.abs(v)) : v;

	}
	public void setMaxDistance(double maxDistance) {
		this.maxDistance = maxDistance;
	}
    
    public void setTemperature(double temp) {
        weightCluster.setTemperature(temp);
    }
}
