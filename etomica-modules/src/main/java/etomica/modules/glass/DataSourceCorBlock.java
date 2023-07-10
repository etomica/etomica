/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.integrator.IntegratorMD;

public class DataSourceCorBlock extends DataSourceBlockAvgCor implements IDataSourceCorBlock {

    protected final double[][] blockSum, corSum, corSelfSum, blockSumIJ;
    protected final long[] nSamples;

    // DataSourceMSD (by default) will collect interval<3, but not every time.  Our superclass
    // doesn't understand that.

    public DataSourceCorBlock(IntegratorMD integrator) {
        super(integrator);
        nSamples = new long[40];
        blockSum = new double[40][40];
        corSum = new double[40][40];
        corSelfSum = new double[40][40];
        blockSumIJ = new double[40][40];
    }

    public void putBlock(int interval, long step, double block) {
        if (interval < minInterval) return;
        processData(interval, step, block);


        nSamples[interval]++;
        for (int k=0; k<blockSum[interval].length; k++) {
            // sum this now for computing correlation with longer blocks
            blockSum[interval][k] += block;
        }
        for (int k=0; k<=interval; k++) {
            // contribution for correlation of i with (shorter block) k
            int y = interval-Math.max(k,minInterval);
            double z = blockSum[k][interval-k]/(1<<y);
            corSum[interval][k] += blockSum[interval][0] * z;
            corSelfSum[interval][k] += z*z;
            blockSumIJ[k][interval-k] += z;
            blockSum[k][interval-k] = 0;
        }
    }

    public double[][] getFullCorrelation() {

        int nData = 0;
        for (int i=minInterval; i<nSamples.length; i++) {
            if (nSamples[i] < 2) {
                nData = i;
                break;
            }
        }
        double[][] cor = new double[nData][nData];
        for (int i = minInterval; i < cor.length; i++) {
            double ai = blockSumIJ[i][0]/nSamples[i];
            double xi = corSelfSum[i][i]/nSamples[i] - ai*ai;
            for (int j = minInterval; j <= i; j++) {
                double aj = blockSumIJ[j][i-j]/nSamples[i];
                double xj = corSelfSum[i][j]/nSamples[i] - aj*aj;
                double xij = corSum[i][j]/nSamples[i] - ai*aj;
                cor[j][i] =  cor[i][j] = xij / Math.sqrt(xi*xj);
            }
        }
        return cor;

    }

}
