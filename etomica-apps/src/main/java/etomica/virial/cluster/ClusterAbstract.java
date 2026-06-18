/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;


import etomica.virial.BoxCluster;

public interface ClusterAbstract {

    /**
     * Returns another instance of an identical cluster (shallow copy). 
     */
    public ClusterAbstract makeCopy();
    
    /**
     * Number of points in the cluster.
     * @return int
     */
    public int pointCount();
    
    /**
     * Value of this cluster for the given pairset at the specified reciprocal
     * temperature.
     */
    public double value(BoxCluster box);

    public void setTemperature(double temperature);
    
}
