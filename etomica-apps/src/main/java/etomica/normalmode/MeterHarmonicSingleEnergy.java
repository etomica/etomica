/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.space.Vector;
import etomica.units.dimensions.Energy;

/**
 * Meter that calculates the Boltzmann-factored harmonic energy of each normal mode for a 
 * configuration given eigenvectors and omegas corresponding to wave vectors.
 * 
 * @author Andrew Schultz
 */
public class MeterHarmonicSingleEnergy implements IDataSource {

    public MeterHarmonicSingleEnergy(CoordinateDefinition coordinateDefinition, NormalModes normalModes) {
        this.coordinateDefinition = coordinateDefinition;
        this.normalModes = normalModes;
        dataInfo = new DataInfoDoubleArray("Harmonic single energy", Energy.DIMENSION, new int[]{0});
        tag = new DataTag();
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    

    public IData getData() {
        double[] x = data.getData();
        
        for (int iVector = 0; iVector < waveVectors.length; iVector++) {
            coordinateDefinition.calcT(waveVectors[iVector], realT, imaginaryT);
            // we want to calculate Q = A T
            // where A is made up of eigenvectors as columns
            int coordinateDim = coordinateDefinition.getCoordinateDim();
            for (int i=0; i<coordinateDim; i++) {
                if (Double.isInfinite(omegaSquared[iVector][i])) {
                    continue;
                }
                double realCoord = 0, imaginaryCoord = 0;
                for (int j=0; j<coordinateDim; j++) {
                    realCoord += realT[j] * eigenvectors[iVector][i][j];
                    imaginaryCoord += imaginaryT[j] * eigenvectors[iVector][i][j];
                }
                double normalCoord = (realCoord*realCoord + imaginaryCoord*imaginaryCoord);
                x[iVector*coordinateDim+i] = waveVectorCoefficients[iVector] * 
                        normalCoord * omegaSquared[iVector][i];
            }
        }
        return data;
    }
    
    public Box getBox() {
        return coordinateDefinition.getBox();
    }

    public void setBox(Box newBox) {
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        
        normalModes.getWaveVectorFactory().makeWaveVectors(newBox);
        setWaveVectors(normalModes.getWaveVectorFactory().getWaveVectors(),normalModes.getWaveVectorFactory().getCoefficients());
        setEigenvectors(normalModes.getEigenvectors());
        setOmegaSquared(normalModes.getOmegaSquared());

        dataInfo = new DataInfoDoubleArray("Harmonic single energy", Energy.DIMENSION, new int[]{waveVectors.length,coordinateDim});
        data = new DataDoubleArray(new int[]{waveVectors.length,coordinateDim});


        realT = new double[coordinateDim];
        imaginaryT = new double[coordinateDim];
    }
    
    public void setWaveVectors(Vector[] newWaveVectors, double[] coefficients) {
        waveVectors = newWaveVectors;
        waveVectorCoefficients = coefficients;
    }
    
    public void setEigenvectors(double[][][] eigenVectors) {
        this.eigenvectors = eigenVectors.clone();
    }
    
    public void setOmegaSquared(double[][] omega2) {
        omegaSquared = new double[omega2.length][omega2[0].length];
        for (int i=0; i<omegaSquared.length; i++) {
            for (int j=0; j<omegaSquared[i].length; j++) {
                // omega is sqrt(kT)/eigenvalue
                omegaSquared[i][j] = omega2[i][j];
            }
        }
    }
    
    public void setName(String newName) {
        name = newName;
    }
    
    public String getName() {
        return name;
    }
    
    private static final long serialVersionUID = 1L;
    protected final CoordinateDefinition coordinateDefinition;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    private final DataTag tag;
    protected double[] realT, imaginaryT;
    protected Vector[] waveVectors;
    protected double[] waveVectorCoefficients;
    protected double[][][] eigenvectors;
    protected double[][] omegaSquared;
    protected String name;
    protected NormalModes normalModes;
}
