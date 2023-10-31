/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;


import etomica.virial.BoxCluster;

/**
 * @author andrew
 *
 * cluster weight wrapper (absolute value of wrapped cluster)
 */
public class ClusterWeightAbs implements ClusterWeight, java.io.Serializable {
	
    protected final ClusterAbstract weightCluster;
    protected boolean doAbs = true;
	
	public ClusterWeightAbs(ClusterAbstract cluster) {
		weightCluster = cluster;
	}

    /**
     * Sets the cluster to actually return the |value| (as default) or just value.
     */
    public void setDoAbs(boolean doAbs) {
        this.doAbs = doAbs;
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
		double v = weightCluster.value(box);
		return doAbs ? Math.abs(v) : v;
	}
    
    public void setTemperature(double temp) {
        weightCluster.setTemperature(temp);
    }
}
