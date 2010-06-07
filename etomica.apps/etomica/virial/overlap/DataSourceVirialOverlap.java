package etomica.virial.overlap;

import etomica.data.AccumulatorRatioAverage;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.util.numerical.AkimaSpline;

/**
 * Measures ratio of two cluster integrals using overlap sampling.  The resulting ratio is 
 * formed from the ratio of a target and reference ratio from two different sub-simulations. 
 */
public class DataSourceVirialOverlap {

    private static final long serialVersionUID = 1L;
    private AccumulatorVirialOverlapSingleAverage refAccumulator, targetAccumulator;
	private final int nBennetPoints;
	
	public DataSourceVirialOverlap(AccumulatorVirialOverlapSingleAverage aRefAccumulator, 
			AccumulatorVirialOverlapSingleAverage aTargetAccumulator) {
		refAccumulator = aRefAccumulator;
		targetAccumulator = aTargetAccumulator;
		nBennetPoints = aRefAccumulator.getNBennetPoints();
	}
    
    public AccumulatorVirialOverlapSingleAverage[] getAccumulators() {
        return new AccumulatorVirialOverlapSingleAverage[]{refAccumulator, targetAccumulator};
    }

    /**
     * Returns the ratio of the reference to target overlap-to-virial ratios
     * (which reduces to reference/target) for the optimal value of the Bennett
     * parameter.  The optimal value is determined via interpolation
     * of y=(ln(ratio) - ln(alpha)) vs. x=ln(alpha) and taking the x value
     * where y=0.
     * 
     * If there are at least 5 alpha values, Akima spline interpolation is
     * used.  With 1 alpha value, no interpolation is needed.  With 3 alpha
     * values, linear interpolation is used.  
     */
	public double[] getOverlapAverageAndError() {
	    if (nBennetPoints == 1) {
	        return new double[]{getAverage(0), getError(0)};
	    }

        double[] lnAlpha = new double[nBennetPoints];
        double[] err = new double[nBennetPoints];
        double[] lnAlphaDiff = new double[nBennetPoints];

        for (int j=0; j<nBennetPoints; j++) {
            lnAlphaDiff[j] += Math.log(getAverage(j));

            double jAlpha = refAccumulator.getBennetBias(j);
            err[j] = getError(j);
            lnAlpha[j] = Math.log(jAlpha);
            lnAlphaDiff[j] -= lnAlpha[j];
        }

        double newAlpha = 0, newErr = 0;
        if (lnAlphaDiff[0] < 0) {
            // first new alpha is less than initial first alpha
            newAlpha = Math.exp(lnAlphaDiff[0] + lnAlpha[0]);
            newErr = newAlpha*err[0];
        }
        else if (lnAlphaDiff[nBennetPoints-1] > 0) {
            newAlpha = Math.exp(lnAlphaDiff[nBennetPoints-1] + lnAlpha[nBennetPoints-1]);
            newErr = newAlpha*err[nBennetPoints-1];
        }
        else if (nBennetPoints > 4) {
            AkimaSpline spline = new AkimaSpline();
            spline.setInputData(lnAlpha, lnAlphaDiff);
            double min = lnAlpha[0];
            double max = lnAlpha[nBennetPoints-1];
            double ymin = lnAlphaDiff[0];
            double ymax = lnAlphaDiff[nBennetPoints-1];
            double[] x = new double[1];
            x[0] = min+(max-min)/(ymax-ymin)*(-ymin);
            while (x[0] > min && x[0] < max) {
                x[0] = 0.5 * (min + max);
                double y = spline.doInterpolation(x)[0];
                if (y == 0) {
                    break;
                }
                if (y*min > 0) {
                    min = x[0];
                    ymin = y;
                }
                else {
                    max = x[0];
                    ymax = y;
                }
                x[0] = min+(max-min)/(ymax-ymin)*(-ymin);
            }
            spline.setInputData(lnAlpha, err);
            newErr = spline.doInterpolation(x)[0];
            newAlpha = Math.exp(x[0]);
        }
        else {
            //linear interpolation (only 3 points)
            for (int i=0; i<nBennetPoints; i++) {
                if (lnAlphaDiff[i] > 0 && lnAlphaDiff[i+1] < 0) {
                    double ix = lnAlpha[i]+(lnAlpha[i+1]-lnAlpha[i])/(lnAlphaDiff[i+1]-lnAlphaDiff[i])*(-lnAlphaDiff[i]);
                    newAlpha = Math.exp(ix);
                    newErr = err[i] + (err[i+1]-err[i])*(ix-lnAlpha[i])/(lnAlpha[i+1]-lnAlpha[i]);
                }
            }
        }

        return new double[]{newAlpha, newErr};
	}
    
    /**
     * Returns the ratio of the reference to target overlap-to-virial ratios
     * (which reduces to target/reference) for the given value of the Bennet
     * parameter.
     */
	public double getAverage(int iParam) {
        double targetAvg = ((DataDoubleArray)((DataGroup)targetAccumulator.getData(iParam)).getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1];
        double refAvg = ((DataDoubleArray)((DataGroup)refAccumulator.getData(iParam)).getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1];
        return refAvg/targetAvg;
	}
	
	/**
     * Returns the index of the Bennet parameter which minimizes the differences 
     * between the Bennet "sums" from the target and references accumulators.  This
     * parameter should be optimal for overlap sampling.
	 */
    public int minDiffLocation() {
		int minDiffLoc = 0;
        double ratio = refAccumulator.getBennetAverage(0)/targetAccumulator.getBennetAverage(0);
        double bias = refAccumulator.getBennetBias(0);
		double minDiff = ratio/bias + bias/ratio - 2;
		for (int i=1; i<nBennetPoints; i++) {
            ratio = refAccumulator.getBennetAverage(i)/targetAccumulator.getBennetAverage(i);
            bias = refAccumulator.getBennetBias(i);
            double newDiff = ratio/bias + bias/ratio - 2;
			if (newDiff < minDiff) {
				minDiffLoc = i;
				minDiff = newDiff;
			}
		}
		return minDiffLoc;
	}
	
    /**
     * Returns the error in the ratio of the reference to target 
     * overlap-to-virial ratios (which reduces to target/reference) 
     * for the optimal value of the Bennet parameter.
     */
	public double getError() {
		return getError(minDiffLocation());
	}

    /**
     * Returns the error in the ratio of the reference to target 
     * overlap-to-virial ratios (which reduces to target/reference) 
     * for the given value of the Bennet parameter.
     */
	public double getError(int iParam) {
		double avg = getAverage(iParam);
        DataGroup dataGroup = (DataGroup)refAccumulator.getData(iParam);
		double refErr = ((DataDoubleArray)dataGroup.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1];
        double refAvg = ((DataDoubleArray)dataGroup.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1];
        double refRelErr = refErr/refAvg;
        dataGroup = (DataGroup)targetAccumulator.getData(iParam);
		double targetErr = ((DataDoubleArray)dataGroup.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1];
		double targetAvg = ((DataDoubleArray)dataGroup.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1];
        double targetRelErr = targetErr/targetAvg;
		return Math.abs(avg)*Math.sqrt(refRelErr*refRelErr+targetRelErr*targetRelErr);
	}

    /**
     * Convenience method that resets reference and target accumulators
     */
    public void reset() {
        targetAccumulator.reset();
        refAccumulator.reset();
    }
}