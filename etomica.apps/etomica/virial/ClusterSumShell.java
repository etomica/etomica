package etomica.virial;

public class ClusterSumShell extends ClusterSum {

    /**
     * Constructs a cluster (shell) that depends on another ClusterSum instance
     * (the core) to do the actual work of calculating all the f values.  This
     * class still does the work of multiplying f bonds together and summing
     * diagrams together.
     * 
     * MayerFunction array is needed only to make the ClusterSum superclass
     * happy.  Length of array should probably match core's f array.
     */
    public ClusterSumShell(ClusterSum coreCluster, ClusterBonds[] subClusters,
            double[] subClusterWeights, MayerFunction[] f) {
        super(subClusters, subClusterWeights, f);
        // grab this now and we won't need to do this again
        fValues = coreCluster.getFValues();
        this.coreCluster = coreCluster;
    }
    
    public ClusterAbstract makeCopy() {
        ClusterSumShell copy = new ClusterSumShell(coreCluster, clusters, clusterWeights, f);
        copy.setTemperature(1/beta);
        return copy;
    }
    
    protected void updateF(BoxCluster box) {
        // our f values are always up to date
    }
    
    protected final ClusterSum coreCluster;
}
