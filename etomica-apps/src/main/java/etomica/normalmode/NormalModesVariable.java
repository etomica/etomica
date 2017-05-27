/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.api.IVector;
import etomica.space.Space;

public class NormalModesVariable {

    public NormalModesVariable(Space space, int nModes, CoordinateDefinition coordinateDefinition) {
        eigenVectors = new double[nModes][coordinateDefinition.getCoordinateDim()];
        waveVectors = new IVector[nModes];
        for (int i=0; i<nModes; i++) {
            waveVectors[i] = space.makeVector();
        }
        phaseAngles = new double[nModes];
    }
    
    public double[][] getEigenVectors() {
        return eigenVectors;
    }

    public IVector[] getWaveVectors() {
        return waveVectors;
    }
    
    public double[] getPhaseAngles() {
        return phaseAngles;
    }

    protected IVector[] waveVectors;
    protected double[][] eigenVectors;
    protected double[] phaseAngles;
}
