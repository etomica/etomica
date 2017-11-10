/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic.overlap;

import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSink;
import etomica.data.types.DataDoubleArray;

/**
 * Accumulator for taking ratio between two sums (and pretend it's an "average")
 * Data added to this accumulator must be a 2-element DataDoubleArray
 */
public class AccumulatorOverlapAverageUa2 implements IDataSink {

    public AccumulatorOverlapAverageUa2(int aNBennetPoints, boolean aIsReference) {
		super();
        isReference = aIsReference;
        setNumAlpha(aNBennetPoints);
        setBennetParam(1.0,5);
        blockCounter = 0;
    }

    public void putDataInfo(IDataInfo info) {
    }
	
	/**
	 * sets the range of parameter values used for Bennets method.
	 * Default is a span of 5 centered about 1 (exp(-5) to (exp(5)).
	 * @param aCenter geometric mean of all values
	 * @param aSpan natural log of ratio of max value to aCenter
	 */
	public void setBennetParam(double aCenter, double aSpan) {
        if (aSpan < 0.0 || (aSpan == 0 && nBennetPoints > 1) || aCenter <= 0.0 ) throw new IllegalArgumentException("span and center must be positive");
        alphaSpan = aSpan;
        alphaCenter = aCenter;
		if (nBennetPoints==1) {
			expX[0] = aCenter;
			return;
		}
		for (int i=0; i<nBennetPoints; i++) {
			expX[i] = Math.exp(2.0*aSpan*(i-(nBennetPoints-1)/2)/(nBennetPoints-1))*aCenter;
		}
		reset();
	}
	
	public void setNumAlpha(int newNumAlpha) {
	    nBennetPoints = newNumAlpha;
        if (nBennetPoints%2 == 0) {
            throw new IllegalArgumentException("try again with an odd aNPoints");
        }
        overlapSum = new double[nBennetPoints];
        overlapSumSquare = new double[nBennetPoints];
        expX = new double[nBennetPoints];
        if (alphaCenter > 0) {
            setBennetParam(alphaCenter, alphaSpan);
        }
        else {
            reset();
        }
	}
	
	public double getAlphaCenter() {
	    return alphaCenter;
	}
	
	public double getAlphaSpan() {
	    return alphaSpan;
	}

    /**
     * Add the given values to the sums and block sums
     */
    // the values coming in should be gamma1/|gamma1| and |gamma2|/|gamma1|
    // where 1 and 2 are target and reference or vica versa, depending on
    // which box the values are coming from.
    public void putData(IData value) {
        double value1 = ((DataDoubleArray)value).getData()[1];
        for (int j=0; j<nBennetPoints; j++) {
            double v;
            // this is actually blockSum[1], but for all the various values of the overlap parameter
            // this doesn't look right, but it is.
            // http://rheneas.eng.buffalo.edu/~andrew/overlapf.pdf
            v = 1;
            if (isReference) {
                v /= (1.0 + expX[j]/value1);
            }
            else {
                v /= (expX[j] + 1.0/value1);
            }
            v = v*v;
            overlapSum[j] += v;
            overlapSumSquare[j] += v*v;
        }
        count++;
    }
    
    public int getNBennetPoints() {
    	return nBennetPoints;
    }
    
    /**
     * Returns iParam'th factor used in the Bennet sum.  Higher values
     * indicate the overlap function is more like the 0th value.
     */
    public double getBennetBias(int iParam) {
        return expX[iParam];
    }
    
    /**
     * Implements DataSource interface, but you probably want to
     * getData for a specific Bennet parameter.
     */
    public IData getData() {
        return null;
    }
    
    /**
     * Returns average value of expression used to determine optimal
     * Bennet parameter (value[1]/(value[1]+expX[iParam]).
     */
    public double getUa2(int iParam) {
        return overlapSum[iParam]/count;
    }

 	/**
 	 * Resets all sums to zero
 	 */
    public void reset() {
        for (int i=0; i<nBennetPoints; i++) {
            overlapSum[i] = 0.0;
            overlapSumSquare[i] = 0.0;
        }
    }
    
    private static final long serialVersionUID = 1L;
    private double[] overlapSumSquare;
    private double[] overlapSum;
    private int nBennetPoints;
    private double[] expX;
    private final boolean isReference;
    protected long blockCounter;
    protected long count;
    protected double alphaCenter, alphaSpan;

}
