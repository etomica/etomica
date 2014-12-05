/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import java.util.Arrays;

/**
 * repulsive potential: Lennard-Jones potential
 * attractive potential: short-handed square well potential 
 * Wertheim's double attraction-site model, hyper-point diagram
 * 2 F-bond, 1 f_R bond diagram
 * @author Hye Min Kim
 */


public class ClusterSumAssociation extends ClusterSum {

    /**
     * Constructor for ClusterSum.
     */
    public ClusterSumAssociation(ClusterBonds[] subClusters, double[] subClusterWeights, MayerFunction[] fArray, ClusterCriteria clusterCriteria) {
        super(subClusters, subClusterWeights, fArray);
        this.clusterCriteria = clusterCriteria;
        
    }

    
    public ClusterAbstract makeCopy() {
        ClusterSumAssociation copy = new ClusterSumAssociation(clusters,clusterWeights,f,clusterCriteria);
        copy.setTemperature(1/beta);
        return copy;
    }

    
    protected void calcValue() {
    	if (clusterCriteria.isValid(fValues)){
    		super.calcValue();  	
        } else {
        	value = 0.0;
        }
    }
    
    public interface ClusterCriteria {
    	public boolean isValid(double [][][]fValues);
    }
    private static final long serialVersionUID = 1L;
    protected final ClusterCriteria clusterCriteria;
}
