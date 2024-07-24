/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.overlap;

import etomica.data.AccumulatorAverage.StatType;
import etomica.data.AccumulatorRatioAverageCovarianceFull;
import etomica.data.IData;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.overlap.AlphaSource;
import etomica.overlap.IntegratorOverlap.ReferenceFracSource;
import etomica.util.Debug;
import etomica.math.numerical.AkimaSpline;

/**
 * Convenience class that analyzes data from the reference and target system,
 * providing average and error estimates on values and ratios.  With multiple
 * alpha values, linear or spline interpolation is used to interpolate and find
 * results for intermediate alpha values.
 * 
 * The targetIndex parameter taken by some methods refers to the target value
 * of interest.  Typically there is one target value, so targetIndex=0.  If
 * more than one target exists (perhaps a nominal target as well as values for
 * specific diagrams, or values for differences of integrands), then those
 * average and ratios can be accessed with the appropriate targetIndex.  The
 * appropriate targetIndex depends on the order that values were returned by
 * the MeterVirial.  The default targetIndex is 0.
 */

/*
 * where this class retrieves data from the target accumulator,
 * 0 to nTargets-1 are target values, nTargets to nTargets+numAlpha-1 are overlap values
 * when the ratio (or error on the ratio) is needed, the ratio of value i to j
 * is at index i*nTotal+j
 */
public class DataVirialOverlap implements ReferenceFracSource {

    protected final AccumulatorRatioAverageCovarianceFull refAccumulator, targetAccumulator;
    protected final AlphaSource alphaSource;
    protected final AkimaSpline spline;
    protected final FullResult fullOverlapResult, fullSystemResult, fullRatioResult;
    protected final StatType avgStat, errorStat, ratioStat, ratioErrorStat;
    protected boolean doIgnoreRefAvg;

	public DataVirialOverlap(AlphaSource alphaSource, AccumulatorRatioAverageCovarianceFull aRefAccumulator, 
			AccumulatorRatioAverageCovarianceFull aTargetAccumulator) {
		refAccumulator = aRefAccumulator;
		targetAccumulator = aTargetAccumulator;
		avgStat = refAccumulator.AVERAGE;
        errorStat = refAccumulator.ERROR;
        ratioStat = refAccumulator.RATIO;
        ratioErrorStat = refAccumulator.RATIO_ERROR;
		this.alphaSource = alphaSource;
		spline =  new AkimaSpline();
		fullOverlapResult = new FullResult();
        fullSystemResult = new FullResult();
        fullRatioResult = new FullResult();
	}

    public AccumulatorRatioAverageCovarianceFull[] getAccumulators() {
        return new AccumulatorRatioAverageCovarianceFull[]{refAccumulator, targetAccumulator};
    }

    public AlphaSource getAlphaSource() {
        return alphaSource;
    }

    protected int getNumTargets() {
        return ((DataInfoGroup)targetAccumulator.getDataInfo()).getSubDataInfo(0).getLength() - alphaSource.getNumAlpha();
    }

    /**
     * Directs this class to ignore the reference average when deciding how
     * much time to spend on each system.  The reference average is generally
     * the same across all simulations meaning that it can be precomputed.
     * This is helpful when the reference average is difficult (compared to the
     * other averages).  What we need to do for that is to not include the
     * error contribution from the reference average in getIdealRefFraction()
     */
    public void setIgnoreReferenceAvg(boolean newDoIgnoreRefAvg) {
        doIgnoreRefAvg = newDoIgnoreRefAvg;
    }

    /**
     * Returns the total result (appropriate ratio of target, reference and
     * overlap averages) and the uncertainty on that value.  Both are
     * calculated for the optimal value of alpha.
     */
    public double[] getAverageAndError() {
        return getAverageAndError(0);
    }

    public double[] getAverageAndError(int targetIndex) {
        double alpha = getOverlapAverage();
        // this populates fullRatioResult
        getFullRatioResultForAlpha(alpha, targetIndex);

        double avg = fullRatioResult.targetAvg/fullRatioResult.refAvg;

        double refErrRatio = fullRatioResult.refErr/fullRatioResult.refAvg;
        double targetErrRatio = fullRatioResult.targetErr/fullRatioResult.targetAvg;
        double errRatio = Math.sqrt(refErrRatio*refErrRatio + targetErrRatio*targetErrRatio);

        return new double[]{avg, Math.abs(avg)*errRatio};
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
    public double getOverlapAverage() {
	    int numAlpha = alphaSource.getNumAlpha();
	    // nTargets is the number of "target" values.  normally, we would just have
	    // the nonimal target, v/|v|.  But we might have additional target values
	    // either for specific diagrams (as when working with flexible molecules)
	    // or when calculating integrand differences (perhaps between temperatures)
	    // the values from the meter will be
	    // <v/|v|>, <v2/|v|>, <v3/|v|>... <o1/|v|>, <o2/|v|>, <o3/|v|>...
	    // where v is the nominal target value, v2, v3... are extra target values
	    // and o1, o2, o3... are overlap values
        int nTargets = ((DataInfoGroup)targetAccumulator.getDataInfo()).getSubDataInfo(0).getLength() - alphaSource.getNumAlpha();
	    if (numAlpha == 1) {
	        return getOverlapAverage(0);
	    }

        double[] lnAlpha = new double[numAlpha];
        double[] lnAlphaDiff = new double[numAlpha];

        for (int j=0; j<numAlpha; j++) {
            double refOverlap = refAccumulator.getData(avgStat).getValue(j+1);
            double targetOverlap = targetAccumulator.getData(avgStat).getValue(j+nTargets);
            lnAlphaDiff[j] += Math.log(refOverlap/targetOverlap);

            double jAlpha = alphaSource.getAlpha(j);
            lnAlpha[j] = Math.log(jAlpha);
            lnAlphaDiff[j] -= lnAlpha[j];
        }

        double newAlpha = 0;
        if (lnAlphaDiff[0] < 0) {
            // first new alpha is less than initial first alpha
            // Math.exp(lnAlphaDiff[0] + lnAlpha[0]) would be the new estimate
            //  for this input alpha, but it's often a poor estimate
            newAlpha = Math.exp(lnAlpha[0]);
        }
        else if (lnAlphaDiff[numAlpha-1] > 0) {
            newAlpha = Math.exp(lnAlpha[numAlpha-1]);
        }
        else if (numAlpha > 4) {
            // interpolate in ln(new)-ln(old) vs. ln(alpha) 
            spline.setInputData(lnAlpha, lnAlphaDiff);
            // we want to find lnDiff = 0
            double min = lnAlpha[0];
            double max = lnAlpha[numAlpha-1];
            double ymin = lnAlphaDiff[0];
            double ymax = lnAlphaDiff[numAlpha-1];
            double[] x = new double[1];
            x[0] = min+(max-min)/(ymax-ymin)*(-ymin);
            while (x[0] > min && x[0] < max) {
                x[0] = 0.5 * (min + max);
                double y = spline.doInterpolation(x)[0];
                if (y == 0) {
                    break;
                }
                if (y*ymin > 0) {
                    min = x[0];
                    ymin = y;
                }
                else {
                    max = x[0];
                    ymax = y;
                }
                x[0] = min+(max-min)/(ymax-ymin)*(-ymin);
            }
            newAlpha = Math.exp(x[0]);
        }
        else {
            //linear interpolation (only 3 points)
            for (int i=0; i<numAlpha; i++) {
                if (lnAlphaDiff[i] > 0 && lnAlphaDiff[i+1] < 0) {
                    double ix = lnAlpha[i]+(lnAlpha[i+1]-lnAlpha[i])/(lnAlphaDiff[i+1]-lnAlphaDiff[i])*(-lnAlphaDiff[i]);
                    newAlpha = Math.exp(ix);
                }
            }
        }

        return newAlpha;
	}

    /**
     * Simply return the reference and target averages and uncertainties.
     */
    // perhaps not so useful, but easy to do.
    public FullResult getFullSystemResult() {
        return getFullSystemResult(0);
    }

    public FullResult getFullSystemResult(int targetIndex) {
        int nTargets = getNumTargets();
        if (targetIndex < 0 || targetIndex > nTargets-1) {
            throw new RuntimeException("invalid target index");
        }
        fullSystemResult.refAvg = refAccumulator.getData(avgStat).getValue(targetIndex);
        fullSystemResult.targetAvg = targetAccumulator.getData(avgStat).getValue(targetIndex);
        fullSystemResult.refErr = refAccumulator.getData(errorStat).getValue(0);
        fullSystemResult.targetErr = targetAccumulator.getData(errorStat).getValue(0);
        return fullSystemResult;
    }

    /**
     * Returns the overlap average and uncertainty for a given value of alpha.
     * This is useful to get an error estimate on a new value of alpha.
     */
    public double[] getOverlapAverageAndErrorForAlpha(double alpha) {
        getFullOverlapResultForAlpha(alpha);
        double avg = fullOverlapResult.refAvg / fullOverlapResult.targetAvg;
        double refErr = fullOverlapResult.refErr / fullOverlapResult.refAvg;
        double targetErr = fullOverlapResult.targetErr / fullOverlapResult.targetAvg;
        double relErr2 = refErr*refErr + targetErr*targetErr;
        double err = Math.sqrt(avg*avg* relErr2);
        return new double[]{avg, err};
    }
    
    /**
     * Returns the full set of overlap results (averages and uncertainties for
     * target and reference overlap averages).  All results are calculated for
     * the given value of alpha.
     */
    public FullResult getFullOverlapResultForAlpha(double alpha) {
        int nTargets = getNumTargets();
	    int numAlpha = alphaSource.getNumAlpha();
	    if (numAlpha == 1) {
	        fullOverlapResult.refAvg = refAccumulator.getData(avgStat).getValue(1);
	        fullOverlapResult.targetAvg = targetAccumulator.getData(avgStat).getValue(nTargets);
	        fullOverlapResult.refErr = refAccumulator.getData(errorStat).getValue(1);
	        fullOverlapResult.targetErr = targetAccumulator.getData(errorStat).getValue(nTargets);
            return fullOverlapResult;
	    }

        double[] lnAlpha = new double[numAlpha];
        double[][] ratios = new double[2][numAlpha];
        double[][] err = new double[2][numAlpha];

        for (int j=0; j<numAlpha; j++) {
            double jAlpha = alphaSource.getAlpha(j);
            ratios[0][j] = refAccumulator.getData(avgStat).getValue(j+1);
            ratios[1][j] = targetAccumulator.getData(avgStat).getValue(j+nTargets);
            err[0][j] = refAccumulator.getData(errorStat).getValue(j+1);
            err[1][j] = targetAccumulator.getData(errorStat).getValue(j+nTargets);
            lnAlpha[j] = Math.log(jAlpha);
        }

        if (numAlpha > 4) {
            double[] x = new double[]{Math.log(alpha)};
            spline.setInputData(lnAlpha, err[0]);
            fullOverlapResult.refErr = spline.doInterpolation(x)[0];
            spline.setInputData(lnAlpha, err[1]);
            fullOverlapResult.targetErr = spline.doInterpolation(x)[0];

            spline.setInputData(lnAlpha, ratios[0]);
            fullOverlapResult.refAvg = spline.doInterpolation(x)[0];
            spline.setInputData(lnAlpha, ratios[1]);
            fullOverlapResult.targetAvg = spline.doInterpolation(x)[0];
        }
        else {
            //linear interpolation (only 3 points)
            double lnIAlpha = Math.log(alpha);
            for (int i=0; i<numAlpha; i++) {
                if (lnAlpha[i] < lnIAlpha && lnAlpha[i+1] > lnIAlpha) {
                    double fac = (lnIAlpha-lnAlpha[i])/(lnAlpha[i+1]-lnAlpha[i]);
                    fullOverlapResult.refErr = err[0][i] + (err[0][i+1]-err[0][i])*fac;
                    fullOverlapResult.targetErr = err[1][i] + (err[1][i+1]-err[1][i])*fac;
                    fullOverlapResult.refAvg = ratios[0][i] + (ratios[0][i+1]-ratios[0][i])*fac;
                    fullOverlapResult.targetAvg = ratios[1][i] + (ratios[1][i+1]-ratios[1][i])*fac;
                    break;
                }
            }
        }

        return fullOverlapResult;
    }

    /**
     * Returns the full set of ratio results (reference and target ratios and
     * errors on those ratios).  All results are calculated for the given value
     * of alpha.
     * 
     * The errors in the ratios are calculated using the covariance (between
     * the numerator and denominator) measured during the simulation.
     */
    public FullResult getFullRatioResultForAlpha(double alpha) {
        return getFullRatioResultForAlpha(alpha, 0);
    }
    
    public FullResult getFullRatioResultForAlpha(double alpha, int targetIndex) {
        int numAlpha = alphaSource.getNumAlpha();
        int nTargets = getNumTargets();
        int nTotal = nTargets + alphaSource.getNumAlpha();
        if (numAlpha == 1) {
            // this gives us ratio of ref to refOverlap
            fullRatioResult.refAvg = refAccumulator.getData(ratioStat).getValue(1);
            // this gives us ratio of targetIndex to targetOverlap
            fullRatioResult.targetAvg = targetAccumulator.getData(ratioStat).getValue(targetIndex*nTotal + nTargets);
            fullRatioResult.refErr = refAccumulator.getData(ratioErrorStat).getValue(1);
            fullRatioResult.targetErr = targetAccumulator.getData(ratioErrorStat).getValue(targetIndex*nTotal + nTargets);
            return fullRatioResult;
        }

        double[] lnAlpha = new double[numAlpha];
        double[][] ratios = new double[2][numAlpha];
        double[][] err = new double[2][numAlpha];

        for (int j=0; j<numAlpha; j++) {
            double jAlpha = alphaSource.getAlpha(j);
            // this gives us ratio of ref to jth refOverlap
            ratios[0][j] = refAccumulator.getData(ratioStat).getValue(j+1);
            // this gives us ratio of targetIndex to jth targetOverlap
            ratios[1][j] = targetAccumulator.getData(ratioStat).getValue(targetIndex*nTotal + (nTargets+j));
            err[0][j] = refAccumulator.getData(ratioErrorStat).getValue(j+1);
            err[1][j] = targetAccumulator.getData(ratioErrorStat).getValue(targetIndex*nTotal + (nTargets+j));
            lnAlpha[j] = Math.log(jAlpha);
        }

        if (numAlpha > 4) {
            double[] x = new double[]{Math.log(alpha)};
            spline.setInputData(lnAlpha, err[0]);
            fullRatioResult.refErr = spline.doInterpolation(x)[0];
            spline.setInputData(lnAlpha, err[1]);
            fullRatioResult.targetErr = spline.doInterpolation(x)[0];

            spline.setInputData(lnAlpha, ratios[0]);
            fullRatioResult.refAvg = spline.doInterpolation(x)[0];
            spline.setInputData(lnAlpha, ratios[1]);
            fullRatioResult.targetAvg = spline.doInterpolation(x)[0];
        }
        else {
            //linear interpolation (only 3 points)
            double lnIAlpha = Math.log(alpha);
            for (int i=0; i<numAlpha; i++) {
                if (lnAlpha[i] < lnIAlpha && lnAlpha[i+1] > lnIAlpha) {
                    double fac = (lnIAlpha-lnAlpha[i])/(lnAlpha[i+1]-lnAlpha[i]);
                    fullRatioResult.refErr = err[0][i] + (err[0][i+1]-err[0][i])*fac;
                    fullRatioResult.targetErr = err[1][i] + (err[1][i+1]-err[1][i])*fac;
                    fullRatioResult.refAvg = ratios[0][i] + (ratios[0][i+1]-ratios[0][i])*fac;
                    fullRatioResult.targetAvg = ratios[1][i] + (ratios[1][i+1]-ratios[1][i])*fac;
                    break;
                }
            }
        }

        return fullRatioResult;
    }

	/**
     * Returns the ratio of the reference to target overlap-to-virial ratios
     * (which reduces to target/reference) for the given value of the Bennet
     * parameter.
     */
	public double getOverlapAverage(int iParam) {
        int nTargets = getNumTargets();
        double targetAvg = targetAccumulator.getData(avgStat).getValue(iParam+nTargets);
        double refAvg = refAccumulator.getData(avgStat).getValue(iParam+1);
        return refAvg/targetAvg;
	}
	
	/**
     * Returns the index of the Bennet parameter which minimizes the differences 
     * between the Bennet "sums" from the target and references accumulators.  This
     * parameter should be optimal for overlap sampling.
	 */
	@Deprecated
    public int minDiffLocation() {
        int nTargets = getNumTargets();
        int numAlpha = alphaSource.getNumAlpha();
		int minDiffLoc = 0;
		IData refAvg = refAccumulator.getData(avgStat);
        IData targetAvg = targetAccumulator.getData(avgStat);
        double ratio = refAvg.getValue(1)/targetAvg.getValue(nTargets);
        double bias = alphaSource.getAlpha(0);
		double minDiff = ratio/bias + bias/ratio - 2;
		for (int i=1; i<numAlpha; i++) {
            ratio = refAvg.getValue(i+1)/targetAvg.getValue(i+nTargets);
            bias = alphaSource.getAlpha(i);
            double newDiff = ratio/bias + bias/ratio - 2;
			if (newDiff < minDiff) {
				minDiffLoc = i;
				minDiff = newDiff;
			}
		}
		return minDiffLoc;
	}

    public double getIdealRefFraction(double oldFrac) {
        double alpha = getOverlapAverage();
        getFullRatioResultForAlpha(alpha);
        
        double refErrorRatio = Math.abs(fullRatioResult.refErr/fullRatioResult.refAvg);
        if (doIgnoreRefAvg) {
            // we've been directed to ignore the reference average.  so, we're just left
            // with the reference overlap average.
            getFullOverlapResultForAlpha(alpha);
            refErrorRatio = Math.abs(fullOverlapResult.refErr/fullOverlapResult.refAvg);
        }
        if (Double.isNaN(refErrorRatio) || refErrorRatio > 1) {
            // if we don't have enough data to calc the error, assume it's 100%
            // if apparent error is > 100%, cap it there (>100% just means our estimate for ratio is bad)
            refErrorRatio = 1.0;
        }

        double targetErrorRatio = Math.abs(fullRatioResult.targetErr/fullRatioResult.targetAvg);
        if (Double.isNaN(targetErrorRatio) || targetErrorRatio > 1) {
            targetErrorRatio = 1.0;
        }
        double refFrac = 1.0 / (1 + targetErrorRatio/refErrorRatio * Math.sqrt((1-oldFrac)/(oldFrac)));
        
        if (Debug.ON && Debug.DEBUG_NOW) {
           // System.out.println("error ratios "+refErrorRatio+" "+targetErrorRatio);
           // System.out.println("frac "+refFrac);
        }
        return refFrac;
    }

    /**
     * Convenience method that resets reference and target accumulators
     */
    public void reset() {
        targetAccumulator.reset();
        refAccumulator.reset();
    }

    public static class FullResult {
        public double refAvg, refErr;
        public double targetAvg, targetErr;
    }
}
