/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;


/**
 * Umbrella cluster wraps a ClusterWheatleySoftDerivatives.  The umbrella cluster's
 * value is the weighted sum of the absolute value of the derivative values (0-n).
 * The weights are equal by default, but can be set.
 */
public class ClusterDerivativeUmbrella implements ClusterWeight, java.io.Serializable {

    /**
     * Contructs an umbrella cluster that wraps a derivatives cluster
     */
    public ClusterDerivativeUmbrella(ClusterWheatleySoftDerivatives cluster) {
        this.cluster = cluster;
        int n = cluster.nDer + 1;
        weightCoefficients = new double[n];
        for (int i = 0; i < n; i++) {
            weightCoefficients[i] = 1.0 / weightCoefficients.length;
        }
    }

    public ClusterAbstract makeCopy() {
        ClusterDerivativeUmbrella newCluster = new ClusterDerivativeUmbrella((ClusterWheatleySoftDerivatives) cluster.makeCopy());
        newCluster.setWeightCoefficients(weightCoefficients);
        return newCluster;
    }

    public int pointCount() {
        return cluster.pointCount();
    }

    public double value(BoxCluster box) {
        cluster.value(box);
        double sum = 0.0;
        int n = cluster.nDer + 1;
        for (int i = 0; i < n; i++) {
            double v = cluster.value[i];
            sum += Math.abs(v) * weightCoefficients[i];
        }
        return sum;
    }

    public void setWeightCoefficients(double[] a) {
        System.arraycopy(a, 0, weightCoefficients, 0, a.length);
    }

    public double[] getWeightCoefficients() {
        return weightCoefficients;
    }

    public void setTemperature(double temp) {
        cluster.setTemperature(temp);
    }

    public ClusterWheatleySoftDerivatives getCluster() {
        return cluster;
    }

    private final ClusterWheatleySoftDerivatives cluster;
    private final double[] weightCoefficients;
}
