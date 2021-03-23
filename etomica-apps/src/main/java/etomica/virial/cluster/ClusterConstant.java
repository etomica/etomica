package etomica.virial.cluster;

import etomica.virial.BoxCluster;

/**
 * This class returns a constant value needed when using random placement of atoms in a box to sample the target system.
 * The constant value is selected such that the reference integral evaluates to 1.
 *
 * @author Arpit Bansal
 */

public class ClusterConstant implements ClusterAbstract{

    @Override
    public ClusterAbstract makeCopy() {
        ClusterConstant cc = new ClusterConstant(points, value);
        return cc;
    }

    @Override
    public int pointCount() {
        return points;
    }

    @Override
    public double value(BoxCluster box) {
        return value;
    }

    @Override
    public void setTemperature(double temperature) {

    }

    public ClusterConstant(int nPoints, double constant){
        points = nPoints;
        value = constant;
    }

    private final double value;
    private final int points;
}
