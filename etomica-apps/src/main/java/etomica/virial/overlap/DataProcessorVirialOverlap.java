/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.overlap;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.overlap.AlphaSource;
import etomica.units.dimensions.Null;

/**
 * DataProcessor which takes in values for v/pi and pi2/pi (from MeterVirial)
 * and calculates the overlap values for various values of alpha.
 * This class can also act as an "AlphaSource"
 */
public class DataProcessorVirialOverlap extends DataProcessor implements AlphaSource {

    public DataProcessorVirialOverlap(int aNumAlpha, boolean aIsReference) {
		super();
        isReference = aIsReference;
        setNumAlpha(aNumAlpha);
        setBennetParam(1.0,5);
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        numIncomingValues = inputDataInfo.getLength();
        if (numIncomingValues < 2) {
            throw new RuntimeException("must have at least 2 values");
        }
        setup();
        return dataInfo;
    }
    
    protected void setup() {
        // recreate our data objects.  new dataInfo will also cause downstream to reset.
        // we need numAlpha elements for the overlap values and 1 for the not-overlap value.
        data = new DataDoubleArray(numAlpha+numIncomingValues-1);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("virial", Null.DIMENSION, new int[]{numAlpha+numIncomingValues-1});
    }

    /**
	 * sets the range of parameter values used for Bennets method.
	 * Default is a span of 5 centered about 1 (exp(-5) to (exp(5)).
	 * @param aCenter geometric mean of all values
	 * @param aSpan natural log of ratio of max value to aCenter
	 */
	public void setBennetParam(double aCenter, double aSpan) {
        if (aSpan < 0.0 || (aSpan == 0 && numAlpha > 1) || aCenter <= 0.0 ) throw new IllegalArgumentException("span and center must be positive");
        alphaSpan = aSpan;
        alphaCenter = aCenter;
		if (numAlpha==1) {
			alpha[0] = aCenter;
			return;
		}
		for (int i=0; i<numAlpha; i++) {
			alpha[i] = Math.exp(2.0*aSpan*(i-(numAlpha-1)/2)/(numAlpha-1))*aCenter;
		}
		// don't really need to recreate these, but we want to notify downstream 
		// passing null as incoming dataInfo because we don't pay attention to it.
		putDataInfo(null);
	}

    public void setNumAlpha(int newNumAlpha) {
        if (newNumAlpha < 1) {
            throw new RuntimeException("# of alpha values must be positive");
        }
        numAlpha = newNumAlpha;
        alpha = new double[numAlpha];
        if (alphaCenter > 0) {
            setBennetParam(alphaCenter, alphaSpan);
            // need to recreate our data objects and notify downstream
            putDataInfo(null);
        }
    }

    public int getNumAlpha() {
        return numAlpha;
    }

    public double getAlphaCenter() {
        return alphaCenter;
    }

    public double getAlphaSpan() {
        return alphaSpan;
    }

    /**
     * take in v1/pi1, pi0/pi1, spit out v1/pi, overlap0, overlap1, overlap2... overlapN
     */
    protected IData processData(IData inputData) {
        double[] x = data.getData();
        for (int i=0; i<numIncomingValues-1; i++) {
            // pass thru the target values
            x[i] = inputData.getValue(i);
        }
        double value1 = inputData.getValue(numIncomingValues-1);
        for (int j=0; j<numAlpha; j++) {
            // this is actually the overlap value for all the various values of the overlap parameter
            // this doesn't look right, but it is.
            // http://rheneas.eng.buffalo.edu/~andrew/overlapf.pdf
            if (isReference) {
                x[j+numIncomingValues-1] = 1.0 / (1.0 + alpha[j]/value1);
            }
            else {
                x[j+numIncomingValues-1] = 1.0 / (alpha[j] + 1.0/value1);
            }
        }
        return data;
    }

    /**
     * Returns ith factor used in the overlap value.
     */
    public double getAlpha(int i) {
        return alpha[i];
    }
    
    protected DataDoubleArray data;
    private int numAlpha;
    private double[] alpha;
    private final boolean isReference;
    protected double alphaCenter, alphaSpan;
    protected int numIncomingValues;
}
