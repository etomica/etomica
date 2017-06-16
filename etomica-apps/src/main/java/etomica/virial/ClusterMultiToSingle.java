package etomica.virial;

public class ClusterMultiToSingle implements ClusterAbstract {

    public ClusterMultiToSingle(ClusterAbstractMultivalue cluster,int m) {
        this.cluster=cluster;
        this.m = m;
               
    }

    @Override
    public ClusterAbstract makeCopy() {        
        return new ClusterMultiToSingle((ClusterAbstractMultivalue)cluster.makeCopy(), m);
    }

    @Override
    public int pointCount() {
        return cluster.pointCount();
    }

    @Override
    public double value(BoxCluster box) {
        return cluster.getAllLastValues(box)[m];
    }

    @Override
    public void setTemperature(double temperature) {
        cluster.setTemperature(temperature);

    }
    protected final ClusterAbstractMultivalue cluster;
    protected final int m;
}
