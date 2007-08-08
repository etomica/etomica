package etomica.virial;


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
