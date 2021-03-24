package etomica.virial;

public class ClusterChainGaussian implements ClusterAbstract{
    protected final double [][] meanPosition;
    protected final double [][] standardDeviation;

    public ClusterChainGaussian (double [][] meanPosition, double [][] standardDeviation) {
        this.meanPosition = meanPosition;
        this.standardDeviation = standardDeviation;
    }
    @Override
    public ClusterAbstract makeCopy() {
        ClusterChainGaussian c = new ClusterChainGaussian(meanPosition, standardDeviation);
        c.setTemperature(1);
        return c;
    }

    @Override
    public int pointCount() {
        return 0;
    }

    @Override
    public double value(BoxCluster box) {
        double value = 1.0;
        for (int i = 0; i < box.getLeafList().size(); i++) {
            for (int j = 0; j < box.getSpace().D(); j++) {
                double a = (box.getLeafList().get(i).getPosition().getX(j) - meanPosition[i][j])/standardDeviation[i][j];
                value *= Math.exp(-a*a/2)/(standardDeviation[i][j]*Math.sqrt(2*Math.PI));
            }
        }
        return value;
    }

    @Override
    public void setTemperature(double temperature) {

    }
}
