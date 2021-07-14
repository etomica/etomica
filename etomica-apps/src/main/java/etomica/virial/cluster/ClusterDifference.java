/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;


import etomica.virial.BoxCluster;

/**
 * This class holds multiple clusters and returns the value from one minus the
 * average value from the others.  This could be implemented as a ClusterSum
 * with a set of ClusterBonds for each of the clusters here, but it's a lot
 * easier to do it this way.
 * 
 * @author Andrew Schultz
 */
public class ClusterDifference implements ClusterAbstract {
    protected final ClusterAbstract clusterAdd;
    protected final ClusterAbstract[] clusterSubtract;

    public ClusterDifference(ClusterAbstract fullTargetCluster, ClusterAbstract[] targetSubtract) {
        this.clusterAdd = fullTargetCluster;
        this.clusterSubtract = targetSubtract;
    }

    public double value(BoxCluster box) {
        // we recalculate everything every time.  we rely on the clusters to
        // cache their values
        double v = clusterAdd.value(box);
        for (int i=0; i<clusterSubtract.length; i++) {
            v -= clusterSubtract[i].value(box)/clusterSubtract.length;
        }
        return v;
    }

    public void setTemperature(double temperature) {
        clusterAdd.setTemperature(temperature);
        for (int i=0; i<clusterSubtract.length; i++) {
            clusterSubtract[i].setTemperature(temperature);
        }
    }

    public int pointCount() {
        return clusterAdd.pointCount();
    }

    public ClusterAbstract makeCopy() {
        ClusterAbstract[] clusterSubtractCopy = new ClusterAbstract[clusterSubtract.length];
        for (int i=0; i<clusterSubtract.length; i++) {
            clusterSubtractCopy[i] = clusterSubtract[i].makeCopy();
        }
        return new ClusterDifference(clusterAdd.makeCopy(), clusterSubtractCopy);
    }
}