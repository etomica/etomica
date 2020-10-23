/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.overlap;

import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.IData;
import etomica.data.types.DataDoubleArray;
import etomica.util.Debug;

import java.io.FileWriter;
import java.io.IOException;

/**
 * Accumulator for taking ratio between two sums (and pretend it's an "average")
 * Data added to this accumulator must be a 2-element DataDoubleArray
 */
public class AccumulatorVirialOverlapSingleAverage extends AccumulatorRatioAverageCovariance {

    public AccumulatorVirialOverlapSingleAverage(int aNBennetPoints, boolean aIsReference) {
		this(1000,aNBennetPoints, aIsReference);
	}
	
    public AccumulatorVirialOverlapSingleAverage(long blockSize, int aNBennetPoints, boolean aIsReference) {
		super(blockSize);
        isReference = aIsReference;
        setNumAlpha(aNBennetPoints);
        setBennetParam(1.0,5);
        blockCounter = 0;
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
        overlapAvg = new double[nBennetPoints];
        blockOverlapAvg = new double[nBennetPoints];
        overlapSumBlockSquare = new double[nBennetPoints];
        overlapSumSquare = new double[nBennetPoints];
        expX = new double[nBennetPoints];
        overlapFirstBlock = new double[nBennetPoints];
        overlapMostRecentBlock = new double[nBennetPoints];
        overlapCorrelationSum = new double[nBennetPoints];
        overlapBlockCovSum = new double[nBennetPoints];
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
    public boolean addData(IData value) {
        if (Debug.ON && value.getLength() != 2) {
            throw new IllegalArgumentException("must receive cluster value and 'other' weight (only)");
        }
        double value1 = ((DataDoubleArray)value).getData()[1];
        for (int j=0; j<nBennetPoints; j++) {
            // this is actually blockSum[1], but for all the various values of the overlap parameter
            // this doesn't look right, but it is.
            // http://rheneas.eng.buffalo.edu/~andrew/overlapf.pdf
            double v;
            if (isReference) {
                v = 1.0 / (1.0 + expX[j]/value1);
            }
            else {
                v = 1.0 / (expX[j] + 1.0/value1);
            }
            blockOverlapAvg[j] += (v - blockOverlapAvg[j]) / (blockSize - blockCountDown + 1);
            overlapSumSquare[j] += v*v;
        }
        // superclass sums up blockSum[1], but we drop it on the floor in doBlockSum in
        // favor of blockOverlapSum
        super.addData(value);
        return true;
    }
    
    protected void doBlockSum() {
    	
    	if (nBennetPoints ==1 && fnm!=null){
    		blockCounter += blockSize;
    		try{
                fileWriter.write(blockCounter + " " + blockOverlapAvg[0] + "\n");
    		
    		} catch (IOException e){
    		
    		}
    	}
    	
        for (int j=0; j<nBennetPoints; j++) {
            double oldAvg = overlapAvg[j];
            overlapAvg[j] += (blockOverlapAvg[j] - overlapAvg[j]) / (count + 1);
			// this is actually blockSum[1], but for all the various values of the overlap parameter
            overlapSumBlockSquare[j] += (blockOverlapAvg[j] - oldAvg) * (blockOverlapAvg[j] - overlapAvg[j]);

            if (!mostRecentBlock.isNaN()) {
                overlapCorrelationSum[j] += overlapMostRecentBlock[j] * blockOverlapAvg[j];
            }
            else {
                overlapFirstBlock[j] = blockOverlapAvg[j];
            }

            overlapBlockCovSum[j] += blockOverlapAvg[j] * currentBlockAvg.getValue(0);

            overlapMostRecentBlock[j] = blockOverlapAvg[j];
            blockOverlapAvg[j] = 0.0;
		}

        super.doBlockSum();
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
        return getData((nBennetPoints-1)/2);
    }
    
    /**
     * Returns average value of expression used to determine optimal
     * Bennet parameter (value[1]/(value[1]+expX[iParam]).
     */
    public double getBennetAverage(int iParam) {
        return overlapAvg[iParam] + (blockOverlapAvg[iParam] - overlapAvg[iParam]) * (blockSize - blockCountDown) / ((double) count * blockSize + (blockSize - blockCountDown));
    }

    /**
     * Return all standard data corresponding to the given Bennet parameter.  
     */
    public IData getData(int iParam) {
        if(count > 0) {
            // fill in data for set "1" with appropriate "overlap" data
            ((DataDoubleArray) average).getData()[1] = overlapAvg[iParam];
            ((DataDoubleArray) currentBlockAvg).getData()[1] = blockOverlapAvg[iParam];
            ((DataDoubleArray) blockVarSum).getData()[1] = overlapSumBlockSquare[iParam];
            ((DataDoubleArray)sumSquare).getData()[1] = overlapSumSquare[iParam];
            ((DataDoubleArray)firstBlock).getData()[1] = overlapFirstBlock[iParam];
            ((DataDoubleArray)mostRecentBlock).getData()[1] = overlapMostRecentBlock[iParam];
            ((DataDoubleArray)correlationSum).getData()[1] = overlapCorrelationSum[iParam];
            double[] x = blockCovSum.getData();
            x[1] = overlapBlockCovSum[iParam];
            x[2] = overlapBlockCovSum[iParam];
            x[3] = overlapSumBlockSquare[iParam];
            // let AccumulatorRatioAverageCovariance do the work for us
            super.getData();
        }
        return dataGroup;
    }
    
 	/**
 	 * Resets all sums to zero
 	 */
    public void reset() {
        for (int i=0; i<nBennetPoints; i++) {
            overlapAvg[i] = 0.0;
            overlapSumBlockSquare[i] = 0.0;
            blockOverlapAvg[i] = 0.0;
            overlapSumSquare[i] = 0.0;
            overlapCorrelationSum[i] = 0.0;
            overlapBlockCovSum[i] = 0;
        }
        super.reset();
    }
    
    public void setFile(String fName){
        this.fnm = fName;
        try {
            fileWriter = new FileWriter(fName);
        } catch (IOException e){
            fileWriter = null;
        }
    }

    public void closeFile(){
        try{
            fileWriter.close();
        } catch (IOException e){}
        fileWriter = null;
    }

    private double[] blockOverlapAvg;
    private double[] overlapSumBlockSquare, overlapSumSquare;
    private double[] overlapAvg;
    protected double[] overlapFirstBlock, overlapMostRecentBlock, overlapCorrelationSum;
    protected double[] overlapBlockCovSum;
    private int nBennetPoints;
    private double[] expX;
    private final boolean isReference;
    protected FileWriter fileWriter;
    protected long blockCounter;
    private String fnm;
    protected double alphaCenter, alphaSpan;

}
