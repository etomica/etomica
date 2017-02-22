/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.space.Space;

/**
 * Obtains wave vectors and coefficients from a file.
 * 
 * @author kofke
 * 
 */
public class WaveVectorFactoryFromFile implements WaveVectorFactory {

    /**
     * @param filename
     *            such that filename.k holds the wave vectors
     * @param D
     *            spatial dimension (not necessarily coordinate dimension)
     */
    public WaveVectorFactoryFromFile(String filename, int D) {
        // read and process wave vectors
        double[][] waveVectorsAndCoefficients = ArrayReader1D
                .getFromFile(filename + ".k");
        waveVectors = new IVectorMutable[waveVectorsAndCoefficients.length];
        coefficients = new double[waveVectors.length];
        double[] justWaveVector = new double[D];
        for (int i = 0; i < waveVectors.length; i++) {
            coefficients[i] = waveVectorsAndCoefficients[i][0];
            for (int j = 0; j < D; j++) {
                justWaveVector[j] = waveVectorsAndCoefficients[i][j + 1];
            }
            Space space = Space.getInstance(justWaveVector.length);
            waveVectors[i] = space.makeVector(justWaveVector);
        }

    }

    public double[] getCoefficients() {
        return coefficients;
    }

    public IVectorMutable[] getWaveVectors() {
        return waveVectors;
    }

    public void makeWaveVectors(IBox box) {
        // nothing to do
    }

    private IVectorMutable[] waveVectors;
    private double[] coefficients;
}
