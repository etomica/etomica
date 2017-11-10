/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.IAction;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

import java.io.Serializable;

/**
 * Calculates the S-matrix for a configuration.  This matrix is formed as T(k) T^(-k), where
 * T is the collective generalized-coordinate vector.
 */
public class MeterNormalMode implements IDataSource, IAction, Serializable {

    public MeterNormalMode() {
        tag = new DataTag();
    }
    
    /**
     * Sets the object that defines the real-space generalized coordinates.
     */
    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
        realT = new double[coordinateDefinition.getCoordinateDim()];
        imaginaryT = new double[coordinateDefinition.getCoordinateDim()];
    }
    
    /**
     * @return the CoordinateDefinition last given via the set method.
     */
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
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
    
    /**
     * Sets the box, and should be called while the Atoms are in 
     * their lattice positions.
     */
    public void setBox(Box newBox) {
        callCount = 0;

        waveVectorFactory.makeWaveVectors(newBox);
        waveVectors = waveVectorFactory.getWaveVectors();
        // we don't actually care about the coefficients
        numWaveVectors = waveVectors.length;

        int coordinateDim = coordinateDefinition.getCoordinateDim();
        DataDoubleArray[] S = new DataDoubleArray[numWaveVectors];
        for (int i=0; i<S.length; i++) {
            // real and imaginary parts
            S[i] = new DataDoubleArray(new int[]{coordinateDim,coordinateDim});
        }
        data = new DataGroup(S);
        DataInfoDoubleArray[] Sinfo = new DataInfoDoubleArray[numWaveVectors];
        CompoundDimension area = new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2});
        for (int i=0; i<Sinfo.length; i++) {
            Sinfo[i] = new DataInfoDoubleArray("S", area, new int[]{coordinateDim,coordinateDim});
        }
        dataInfo = new DataInfoGroup("all S", Null.DIMENSION, Sinfo);
    }
    
    public Box getBox() {
        return coordinateDefinition.getBox();
    }
    
    public Vector[] getWaveVectors() {
        return waveVectors;
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    /**
     * Calculating things and adds terms to the sums
     */
    public void actionPerformed() {
        callCount++;
        int coordinateDim = coordinateDefinition.getCoordinateDim();

        // |data.E(0)| here to calculate the current value rather than the sum
        // loop over wave vectors
        for (int iVector = 0; iVector < numWaveVectors; iVector++) {

            coordinateDefinition.calcT(waveVectors[iVector], realT, imaginaryT);
            
            // add to S(k).  imaginary part of S is 0
            double[] sValues = ((DataDoubleArray)data.getData(iVector)).getData();
            for (int i=0; i<coordinateDim; i++) {
                for (int j=0; j<coordinateDim; j++) {
                    sValues[i*coordinateDim+j] += realT[i]*realT[j] + imaginaryT[i]*imaginaryT[j];
                }
            }
        }
    }

    /**
     * Returns the DataGroup of S(k) Tensors corresponding to the sum of 
     * T(k)*transpose(T(-k)).  To get the average (U), divide by callCount().
     */
    public IData getData() {
        return data;
    }

    /**
     * Sets the tensor summation to 0.
     */
    public void reset() {
        data.E(0);
        callCount = 0;
    }
    
    public int getCallCount() {
        return callCount;
    }
    
    public void setName(String newName) {
        name = newName;
    }
    
    public String getName() {
        return name;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    private static final long serialVersionUID = 1L;
    private Vector[] waveVectors;
    private WaveVectorFactory waveVectorFactory;
    protected CoordinateDefinition coordinateDefinition;
    private int numWaveVectors;
    private String name;
    private final DataTag tag;
    private IDataInfo dataInfo;
    private DataGroup data;
    private int callCount;

    protected double[] realT, imaginaryT;
}
