/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.space.Vector;
import etomica.normalmode.CoordinateDefinition.BasisCell;

/**
 * Class that calculates the dq/dx Jacobian.
 * @author Andrew Schultz
 */
public class CalcJacobian {

    public CalcJacobian() {
    }

    public double[][] getJacobian() {
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        int basisSize = cells[0].molecules.size();
        int l = coordinateDim * cells.length;
        double[][] jacobian = new double[l][l];
        // # of spatial dimensions
        int spaceDim = waveVectors[0].getD();
        
        int vectorPos = 0;
        double sqrtN = Math.sqrt(cells.length);
        for (int iVector = 0; iVector < waveVectors.length; iVector++) {
            double boxAngle = Double.NaN;
            for (int iCell = 0; iCell < cells.length; iCell++) {
                Vector latticePosition = cells[iCell].cellPosition;
                double kR = waveVectors[iVector].dot(latticePosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                if (Math.abs(coskR) < 1.e-14) {
                    coskR = 0;
                }
                if (Math.abs(sinkR) < 1.e-14) {
                    sinkR = 0;
                }
                for (int iDim = 0; iDim < spaceDim; iDim++) {
                    if (waveVectorCoefficients[iVector] == 1) {
                        for (int i=0; i<basisSize; i++) {
                            jacobian[vectorPos*coordinateDim+i*spaceDim+iDim][iCell*coordinateDim+i*spaceDim+iDim] = coskR / sqrtN;
                            jacobian[(vectorPos+1)*coordinateDim+i*spaceDim+iDim][iCell*coordinateDim+i*spaceDim+iDim] = -sinkR / sqrtN;
                        }
                    }
                    else {
                        // single degree of freedom.
                        // All kR must be 100% in-box or 100% out of box. kR-boxAngle will be 0 or pi
                        if (!Double.isNaN(boxAngle) && Math.abs(Math.sin(kR-boxAngle)) > 1.e-10) {
                            throw new RuntimeException("oops "+iVector+" "+waveVectors[iVector]+" "+boxAngle+" "+kR);
                        }
                        boxAngle = kR;
                        double value;
                        // either one works so long as it's not 0
                        if (Math.abs(coskR) > 0.1) {
                            value = coskR > 0 ? 1.0/sqrtN : -1.0/sqrtN;
                        }
                        else {
                            value = sinkR > 0 ? 1.0/sqrtN : -1.0/sqrtN;
                        }
                        for (int i=0; i<basisSize; i++) {
                            jacobian[vectorPos*coordinateDim+i*spaceDim+iDim][iCell*coordinateDim+i*spaceDim+iDim] = value;
                        }
                    }
                }
                // handle orientation here
            }
            vectorPos += Math.round(waveVectorCoefficients[iVector]*2);
        }
        return jacobian;
    }

    /**
     * Sets the object that defines the normal-coordinate wave vectors.
     */
    public void setWaveVectorFactory(WaveVectorFactory newWaveVectorFactory) {
        waveVectorFactory = newWaveVectorFactory;
    }

    /**
     * @return the WaveVectorFactory last given via the set methods.
     */
    public WaveVectorFactory getWaveVectorFactory() {
        return waveVectorFactory;
    }
    
    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
        coordinateDim = coordinateDefinition.getCoordinateDim();

        waveVectorFactory.makeWaveVectors(newCoordinateDefinition.getBox());
        waveVectors = waveVectorFactory.getWaveVectors();
        waveVectorCoefficients = waveVectorFactory.getCoefficients();
    }
    
    private static final long serialVersionUID = 1L;
    protected int coordinateDim;
    protected CoordinateDefinition coordinateDefinition;
    protected Vector[] waveVectors;
    protected double[] waveVectorCoefficients;
    private WaveVectorFactory waveVectorFactory;
}
