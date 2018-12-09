/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.elasticconstantmappedavg;

import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.integrator.mcmove.MCMoveTracker;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveHarmonicmodified {

    public MCMoveHarmonicmodified() { }


    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
        uOld = null;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public void setOmegaSquared(double[][] omega2) {
        oneOverOmega2 = new double[omega2.length][omega2[0].length];
        for (int i=0; i<oneOverOmega2.length; i++) {
            for (int j=0; j<oneOverOmega2[i].length; j++) {
                oneOverOmega2[i][j] = 1.0/Math.sqrt(omega2[i][j]);
            }
        }
    }

    public void setWaveVectorCoefficients(double[] newWaveVectorCoefficients) {
        sqrtWVC = new double[newWaveVectorCoefficients.length];
        for (int i =0; i < newWaveVectorCoefficients.length; i++){
            sqrtWVC[i] = Math.sqrt(2*newWaveVectorCoefficients[i]);
        }
    }

    public void setEigenVectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }

    public void setWaveVectors(Vector[] newWaveVectors) {
        waveVectors = newWaveVectors;
    }

    public double doTrial(double[][] rRand,BasisCell[] cells, int[] modeNum,double[][] omega2) {

        int coordinateDim = coordinateDefinition.getCoordinateDim();
        u = new double[coordinateDim];
        double normalization = 1/Math.sqrt(cells.length);
        double energy=0.0;
        for (int iCell = 0; iCell<cells.length; iCell++) {

            BasisCell cell = cells[iCell];
            for (int i=0; i<coordinateDim; i++) {
                u[i] = 0;
            }
            //loop over wavevectors and sum contribution of each to the generalized coordinates
            for (int iVector=0; iVector<waveVectors.length; iVector++) {
                double kR = waveVectors[iVector].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                     for (int j = 0; j < coordinateDim; j++) {
                        u[j] += eigenVectors[iVector][0][j] * (rRand[iVector][0] * coskR);

                     }
                for (int j = 0; j < coordinateDim; j++) {
if (omega2[iVector][j]==Double.POSITIVE_INFINITY) continue;
                    energy += rRand[iVector][j] * rRand[iVector][j] * omega2[iVector][j];
                }

            }
            for (int i=0; i<coordinateDim; i++) {
                u[i] *= normalization;
            }
            coordinateDefinition.setToU(cell.molecules, u);
        }
return energy;
    }


    protected CoordinateDefinition coordinateDefinition;
    private double[][] oneOverOmega2;
    private double[][][] eigenVectors;
    private Vector[] waveVectors;
    private double[] sqrtWVC;
    protected double[] u;
    protected double[][] rRand;
     protected double temperature;
     protected double[][] uOld;


}
