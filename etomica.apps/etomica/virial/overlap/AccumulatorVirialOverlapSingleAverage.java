package etomica.virial.overlap;

import etomica.data.AccumulatorRatioAverage;

/**
 * Accumulator for taking ratio between two sums (and pretend it's an "average")
 */
public class AccumulatorVirialOverlapSingleAverage extends AccumulatorRatioAverage {

	public AccumulatorVirialOverlapSingleAverage(int aNBennetPoints) {
		super();
		nBennetPoints = aNBennetPoints;
		if (nBennetPoints%2 == 0) {
			throw new IllegalArgumentException("try again with an odd aNPoints");
		}
		bennetSum = new double[nBennetPoints];
		overlapSum = new double[nBennetPoints];
		blockOverlapSum = new double[nBennetPoints];
        overlapSumSquare = new double[nBennetPoints];
        blockOverlapSumSq = new double[nBennetPoints];
		expX = new double[nBennetPoints];
        setBennetParam(1.0,5);
		setNData(2);
	}
	
	/**
	 * sets the range of parameter values used for Bennets method.
	 * Default is a span of 5 centered about 1 (exp(-5) to (exp(5)).  The
     * center for the reference and target cluster should be inverses.
	 * @param aCenter geometric mean of all values (measuring prefernce for 
	 * this cluster compared to the "other" one) 
	 * @param aSpan natural log of ratio of max or min values to aCenter
	 */
	public void setBennetParam(double aCenter, double aSpan) {
		if (nBennetPoints==1) {
			expX[0] = aCenter;
			return;
		}
		if (aSpan <= 0.0) throw new IllegalArgumentException("aSpan must be positive");
		for (int i=0; i<nBennetPoints; i++) {
			expX[i] = Math.exp(2.0*aSpan*(double)(i-(nBennetPoints-1)/2)/(double)(nBennetPoints-1))*aCenter;
            System.out.println("bennet param "+i+" "+expX[i]);
		}
	}

    /**
     * Add the given values to the sums and block sums
     */
    public void addData(double[] value) {
        if(value.length != nData) {
            throw new IllegalArgumentException("must receive cluster value, weight and 'other' weight (only)");
        }
        if (Math.abs(value[0])-1.0 > 1.e-10) {
            throw new IllegalArgumentException("value[0] should be 1 (gamma/weight), but it was "+value[0]);
        }
        for (int j=0; j<nBennetPoints; j++) {
            if (nBennetPoints > 1) {
                bennetSum[j] += value[1]/(value[1]+ expX[j]);
            }
            // this is actually blockSum[1], but for all the various values of the overlap parameter
            double v = value[1]*(1+expX[j]) / (value[1]+expX[j]);
            blockOverlapSum[j] += v;
            blockOverlapSumSq[j] += v*v;
        }
        // superclass sums up blockSum[1], but we drop it on the floor in doBlockSum in
        // favor of blockOverlapSum
        super.addData(value);
    }
    
    protected void doBlockSum() {
        for (int j=0; j<nBennetPoints; j++) {
			// this is actually blockSum[1], but for all the various values of the overlap parameter
            blockOverlapSum[j] /= blockSize;
	    	overlapSum[j] += blockOverlapSum[j];
            overlapSumSquare[j] += blockOverlapSum[j]*blockOverlapSum[j];
            blockOverlapSum[j] = 0.0;
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
    public double[] getData() {
        return getData((nBennetPoints-1)/2);
    }
    
    /**
     * Returns average value of expression used to determine optimal
     * Bennet parameter (value[1]/(value[1]+expX[iParam]).
     */
    public double getBennetAverage(int iParam) {
        return bennetSum[iParam]/((double)count*blockSize+(blockSize-blockCountDown));
    }

    /**
     * Return all standard data corresponding to the given Bennet parameter.  
     */
    public double[] getData(int iParam) {
        int currentBlockCount = blockSize - blockCountDown;
        if(count+currentBlockCount == 0) {
            setNaN(data);
        } else {
            // fill in data for set "1" with appropriate "overlap" data
            sum[1] = overlapSum[iParam];
            blockSum[1] = blockOverlapSum[iParam];
            sumSquare[1] = overlapSumSquare[iParam];
            blockSumSq[1] = blockOverlapSumSq[iParam];
            // let AccumulatorRatioAverage do the work for us
            super.getData();
        }
        return data;
    }
    
 	/**
 	 * Resets all sums to zero
 	 */
     public void reset() {
         for (int i=0; i<nBennetPoints; i++) {
         	bennetSum[i] = 0.0;
         	overlapSum[i] = 0.0;
         	blockOverlapSum[i] = 0.0;
         }
     	 super.reset();
     }
     
    private final double[] overlapSum, blockOverlapSum;
    private final double[] overlapSumSquare, blockOverlapSumSq;
    private final double[] bennetSum;
    private final int nBennetPoints;
    private final double[] expX;
}
