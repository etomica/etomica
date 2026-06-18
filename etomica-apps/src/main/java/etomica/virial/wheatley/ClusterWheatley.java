/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.wheatley;

import etomica.virial.BoxCluster;
import etomica.virial.cluster.ClusterAbstract;

public interface ClusterWheatley extends ClusterAbstract {

    /**
     * Returns edgeCount (number of overlaps) of configuration passed to
     * checkConfig
     */
    public int getEdgeCount();

    /**
     * Returns number of cliques for the configuration passed to
     * checkConfig
     */
    public int getCliqueCount();

    /**
     * Returns number of e-bond cliques for the configuration passed to
     * checkConfig
     */
    public int getECliqueCount();

    /**
     * Returns number bond mask for all points for the configuration passed to
     * checkConfig
     */
    public int[] getFullBondMask();

    /**
     * Returns an array containing the clique sets for the configuration passed
     * to checkConfig.  The array is filled from the beginning and may contain
     * cruft at the end (use getCliqueCount to know where to stop).
     */
    public int[] getCliques();

    /**
     * Returns an array containing the clique sets (based on e-bonds) for the
     * configuration passed to checkConfig.  The array is filled from the
     * beginning and may contain cruft at the end (use getCliqueCount to know
     * where to stop).
     */
    public int[] getECliques();
    
    /**
     * Returns outDegee (number of bonds for each point) of the configuration
     * passed to checkConfig
     */
    public byte[] getOutDegree();

    /**
     * Determines whether the configuration is zero-valued or not.  If not
     * zero, this method also computes various quantities returnable through
     * other methods.
     */
    public boolean checkConfig(BoxCluster box);

    /**
     * Calculates the value of the configuration.  This should be called only
     * after checkConfig.
     */
    public double calcValue(BoxCluster box);
}