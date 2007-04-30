package etomica.virial.overlap;

import etomica.data.AccumulatorRatioAverage;
import etomica.data.Data;
import etomica.data.types.DataDoubleArray;
import etomica.simulation.Simulation;
import etomica.util.Debug;

/**
 * Accumulator for taking ratio between two sums (and pretend it's an "average")
 * Data added to this accumulator must be a 2-element DataDoubleArray
 */
public class AccumulatorVirialOverlapSingleAverage extends AccumulatorRatioAverage {

    public AccumulatorVirialOverlapSingleAverage(Simulation sim, int aNBennetPoints, boolean aIsReference) {
		this(sim.getDefaults().blockSize,aNBennetPoints, aIsReference);
	}
	
    public AccumulatorVirialOverlapSingleAverage(int blockSize, int aNBennetPoints, boolean aIsReference) {
		super(blockSize);
		nBennetPoints = aNBennetPoints;
        isReference = aIsReference;
		if (nBennetPoints%2 == 0) {
			throw new IllegalArgumentException("try again with an odd aNPoints");
		}
		bennetSum = new double[nBennetPoints];
		overlapSum = new double[nBennetPoints];
		blockOverlapSum = new double[nBennetPoints];
        overlapSumSquare = new double[nBennetPoints];
        blockOverlapSumSq = new double[nBennetPoints];
        overlapSumSquareBlock = new double[nBennetPoints];
		expX = new double[nBennetPoints];
        setBennetParam(1.0,5);
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
			expX[i] = Math.exp(2.0*aSpan*(i-(nBennetPoints-1)/2)/(nBennetPoints-1))*aCenter;
//            System.out.println("bennet param "+i+" "+expX[i]);
		}
	}

    /**
     * Add the given values to the sums and block sums
     */
    // the values coming in should be gamma1/|gamma1| and |gamma2|/|gamma1|
    // where 1 and 2 are target and reference or vica versa, depending on
    // which phase the values are coming from.
    public void addData(Data value) {
        if (Debug.ON && ((DataDoubleArray)value).getLength() != 2) {
            throw new IllegalArgumentException("must receive cluster value and 'other' weight (only)");
        }
        double value1 = ((DataDoubleArray)value).getData()[1];
        for (int j=0; j<nBennetPoints; j++) {
            double v;
            if (Double.isInfinite(value1)) {
                if (isReference) {
                    v = 1+expX[j];
                }
                else {
                    v = 1+1.0/expX[j];
                }
            }
            else {
                // this is actually blockSum[1], but for all the various values of the overlap parameter
                // this doesn't look right, but it is.
                // http://rheneas.eng.buffalo.edu/~andrew/overlapf.pdf
                v = 1 + expX[j];
                if (isReference) {
                    v /= (1.0 + expX[j]/value1);
                }
                else {
                    v /= (expX[j] + 1.0/value1);
                }
            }
            bennetSum[j] += v;
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
            overlapSumSquareBlock[j] += blockOverlapSumSq[j];
            blockOverlapSumSq[j] = 0.0;
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
    public Data getData() {
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
    public Data getData(int iParam) {
        if(count > 0) {
            // fill in data for set "1" with appropriate "overlap" data
            ((DataDoubleArray)sum).getData()[1] = overlapSum[iParam];
            ((DataDoubleArray)blockSum).getData()[1] = blockOverlapSum[iParam];
            ((DataDoubleArray)sumSquare).getData()[1] = overlapSumSquare[iParam];
            ((DataDoubleArray)sumSquareBlock).getData()[1] = overlapSumSquareBlock[iParam];
            // let AccumulatorRatioAverage do the work for us
            super.getData();
        }
        return dataGroup;
    }
    
 	/**
 	 * Resets all sums to zero
 	 */
    public void reset() {
        for (int i=0; i<nBennetPoints; i++) {
            bennetSum[i] = 0.0;
            overlapSum[i] = 0.0;
            overlapSumSquare[i] = 0.0;
            blockOverlapSum[i] = 0.0;
            blockOverlapSumSq[i] = 0.0;
            overlapSumSquareBlock[i] = 0.0;
        }
        super.reset();
    }
     
    private static final long serialVersionUID = 1L;
    private final double[] overlapSum, blockOverlapSum;
    private final double[] overlapSumSquare, blockOverlapSumSq;
    private final double[] overlapSumSquareBlock;
    private final double[] bennetSum;
    private final int nBennetPoints;
    private final double[] expX;
    private final boolean isReference;
}
