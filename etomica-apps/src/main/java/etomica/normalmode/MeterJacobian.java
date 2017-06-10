/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import Jama.Matrix;
import etomica.space.Vector;
import etomica.data.DataSourceScalar;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.units.dimensions.Null;

public class MeterJacobian extends DataSourceScalar {

    public MeterJacobian(CoordinateDefinition coordinateDefinition, NormalModesVariable normalModes) {
        super("Jacobian", Null.DIMENSION);
        this.coordinateDefinition = coordinateDefinition;
        this.normalModes = normalModes;
    }

    public double getDataAsScalar() {
        Vector[] waveVectors = normalModes.getWaveVectors();
        double[][] eigenVectors = normalModes.getEigenVectors();
        double[] phaseAngles = normalModes.getPhaseAngles();

        
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        int l = coordinateDim * cells.length;
        double[][] jacobian = new double[waveVectors.length][l];
        // # of spatial dimensions
        
        int vectorPos = 0;
        for (int iMode = 0; iMode < waveVectors.length; iMode++) {
            for (int iCell = 0; iCell < cells.length; iCell++) {
                Vector latticePosition = cells[iCell].cellPosition;
                double kR = waveVectors[iMode].dot(latticePosition);
                double coskR = Math.cos(kR+phaseAngles[iMode]);
                for (int iDim = 0; iDim < coordinateDim; iDim++) {
                    jacobian[iMode][iCell*coordinateDim+iDim] = eigenVectors[iMode][iDim]*coskR/Math.sqrt(cells.length)*Math.sqrt(2);
                }
                // handle orientation here
            }
            vectorPos += 1;
        }
        Matrix m = new Matrix(jacobian);
//        m.print(10, 5);
        return Math.abs(m.det());
    }

    protected final CoordinateDefinition coordinateDefinition;
    protected final NormalModesVariable normalModes;
}
