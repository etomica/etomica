/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.overlap;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.overlap.IntegratorOverlap.ReferenceFracSource;
import etomica.util.Debug;
import etomica.math.numerical.AkimaSpline;

/**
 * Measures ratio of two cluster integrals using overlap sampling.  The resulting ratio is 
 * formed from the ratio of a target and reference ratio from two different sub-simulations. 
 */
public class DataSourceVirialOverlap implements ReferenceFracSource {

    private static final long serialVersionUID = 1L;
    protected final AccumulatorVirialOverlapSingleAverage refAccumulator, targetAccumulator;
	
	public DataSourceVirialOverlap(AccumulatorVirialOverlapSingleAverage aRefAccumulator, 
			AccumulatorVirialOverlapSingleAverage aTargetAccumulator) {
		refAccumulator = aRefAccumulator;
		targetAccumulator = aTargetAccumulator;
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
	    int nBennetPoints = refAccumulator.getNBennetPoints();
	    if (nBennetPoints == 1) {
	        return new double[]{getAverage(0), getError(0)};
	    }

        double[] lnAlpha = new double[nBennetPoints];
        double[] err = new double[nBennetPoints];
        double[] lnAlphaDiff = new double[nBennetPoints];

        for (int j=0; j<nBennetPoints; j++) {
            double refOverlap = ((DataDoubleArray)((DataGroup)refAccumulator.getData(j)).getData(refAccumulator.AVERAGE.index)).getData()[1];
            double targetOverlap = ((DataDoubleArray)((DataGroup)targetAccumulator.getData(j)).getData(refAccumulator.AVERAGE.index)).getData()[1];
            lnAlphaDiff[j] += Math.log(refOverlap/targetOverlap);

            double jAlpha = refAccumulator.getBennetBias(j);
            err[j] = getError(j);
            lnAlpha[j] = Math.log(jAlpha);
            lnAlphaDiff[j] -= lnAlpha[j];
        }

        double newAlpha = 0, newErr = 0;
        if (lnAlphaDiff[0] < 0) {
            // first new alpha is less than initial first alpha
            newAlpha = Math.exp(lnAlpha[0]); //Math.exp(lnAlphaDiff[0] + lnAlpha[0]);
            newErr = Double.NaN; //newAlpha*err[0];
        }
        else if (lnAlphaDiff[nBennetPoints-1] > 0) {
            newAlpha = Math.exp(lnAlpha[nBennetPoints-1]); //Math.exp(lnAlphaDiff[nBennetPoints-1] + lnAlpha[nBennetPoints-1]);
            newErr = Double.NaN; //newAlpha*err[nBennetPoints-1];
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

        double targetAvg = ((DataGroup)targetAccumulator.getData(0)).getData(targetAccumulator.AVERAGE.index).getValue(0); 
        double refAvg = ((DataGroup)refAccumulator.getData(0)).getData(refAccumulator.AVERAGE.index).getValue(0);
        double ratio = targetAvg/refAvg;
        newAlpha *= ratio;
        double refErr = ((DataGroup)refAccumulator.getData(0)).getData(refAccumulator.ERROR.index).getValue(0);
        double targetErr = ((DataGroup)targetAccumulator.getData(0)).getData(targetAccumulator.ERROR.index).getValue(0);
        double refErrRatio = refErr/refAvg;
        double targetErrRatio = targetErr/targetAvg;
        newErr = Math.sqrt(newErr*newErr + ratio*ratio*(refErrRatio*refErrRatio + targetErrRatio*targetErrRatio));

        return new double[]{newAlpha, newErr};
	}

    public double getOverlapAverageForAlpha(double alpha) {
        int nBennetPoints = refAccumulator.getNBennetPoints();

        if (nBennetPoints == 1) {
            throw new RuntimeException("need more than one alpha");
        }

        double[] lnAlpha = new double[nBennetPoints];
        double[] lnRatio = new double[nBennetPoints];

        for (int j=0; j<nBennetPoints; j++) {
            double refOverlap = ((DataDoubleArray)((DataGroup)refAccumulator.getData(j)).getData(refAccumulator.AVERAGE.index)).getData()[1];
            double targetOverlap = ((DataDoubleArray)((DataGroup)targetAccumulator.getData(j)).getData(targetAccumulator.AVERAGE.index)).getData()[1];
            lnRatio[j] = Math.log(refOverlap/targetOverlap);

            double jAlpha = refAccumulator.getBennetBias(j);
            lnAlpha[j] = Math.log(jAlpha);
        }

        double newRatio = 0;
        double lnNewAlpha = Math.log(alpha);
        if (lnNewAlpha < lnAlpha[0] || lnNewAlpha > lnAlpha[nBennetPoints-1]) {
            throw new RuntimeException("alpha value must be within bounds of existing alpha range");
        }
        if (nBennetPoints > 4) {
            AkimaSpline spline = new AkimaSpline();
            spline.setInputData(lnAlpha, lnRatio);
            newRatio = Math.exp(spline.doInterpolation(new double[]{lnNewAlpha})[0]);
        }
        else {
            //linear interpolation (only 3 points)
            for (int i=0; i<nBennetPoints; i++) {
                if (lnAlpha[i] >= lnNewAlpha && lnAlpha[i+1] <= lnNewAlpha) {
                    double ix = lnRatio[i] + (lnRatio[i+1] - lnRatio[i]) / (lnAlpha[i+1] - lnAlpha[i]) * (lnNewAlpha - lnAlpha[i]);
                    newRatio = Math.exp(ix);
                }
            }
        }
        return newRatio;
    }

	/**
     * Returns the ratio of the reference to target overlap-to-virial ratios
     * (which reduces to target/reference) for the given value of the Bennet
     * parameter.
     */
	public double getAverage(int iParam) {
        double targetAvg = ((DataDoubleArray)((DataGroup)targetAccumulator.getData(iParam)).getData(targetAccumulator.RATIO.index)).getData()[1];
        double refAvg = ((DataDoubleArray)((DataGroup)refAccumulator.getData(iParam)).getData(refAccumulator.RATIO.index)).getData()[1];
        return refAvg/targetAvg;
	}
	
	/**
     * Returns the index of the Bennet parameter which minimizes the differences 
     * between the Bennet "sums" from the target and references accumulators.  This
     * parameter should be optimal for overlap sampling.
	 */
    public int minDiffLocation() {
        int nBennetPoints = refAccumulator.getNBennetPoints();
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

    public double getIdealRefFraction(double oldFrac) {
        int minDiffLoc = minDiffLocation();

        DataGroup refData = (DataGroup)refAccumulator.getData(minDiffLoc);
        double refError = ((DataDoubleArray)refData.getData(refAccumulator.RATIO_ERROR.index)).getData()[1];
        double refErrorRatio = refError/Math.abs(((DataDoubleArray)refData.getData(refAccumulator.RATIO.index)).getData()[1]);
        if (Debug.ON && Debug.DEBUG_NOW) {
            System.out.println("0 "+Math.abs(refError)+" "+Math.abs(refError/refErrorRatio));
        }
        if (Double.isNaN(refErrorRatio) || refErrorRatio > 1) {
            // if we don't have enough data to calc the error, assume it's 100%
            // if apparent error is > 100%, cap it there (>100% just means our estimate for ratio is bad)
            refErrorRatio = 1.0;
        }

        DataGroup targetData = (DataGroup)targetAccumulator.getData(minDiffLoc);
        double targetError = ((DataDoubleArray)targetData.getData(targetAccumulator.RATIO_ERROR.index)).getData()[1];
        double targetErrorRatio = targetError/Math.abs(((DataDoubleArray)targetData.getData(targetAccumulator.RATIO.index)).getData()[1]);
        if (Debug.ON && Debug.DEBUG_NOW) {
            System.out.println("1 "+Math.abs(targetError)+" "+Math.abs(targetError/targetErrorRatio));
        }
        if (Double.isNaN(targetErrorRatio) || targetErrorRatio > 1) {
            targetErrorRatio = 1.0;
        }
        double refFrac = 1.0 / (1 + targetErrorRatio/refErrorRatio * Math.sqrt((1-oldFrac)/(oldFrac)));

        if (Debug.ON && Debug.DEBUG_NOW) {
            System.out.println("error ratios "+refErrorRatio+" "+targetErrorRatio);
            System.out.println("frac "+refFrac);
        }
        return refFrac;
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
		double refErr = ((DataDoubleArray)dataGroup.getData(refAccumulator.RATIO_ERROR.index)).getData()[1];
        double refAvg = ((DataDoubleArray)dataGroup.getData(refAccumulator.RATIO.index)).getData()[1];
        double refRelErr = refErr/refAvg;
        dataGroup = (DataGroup)targetAccumulator.getData(iParam);
		double targetErr = ((DataDoubleArray)dataGroup.getData(targetAccumulator.RATIO_ERROR.index)).getData()[1];
		double targetAvg = ((DataDoubleArray)dataGroup.getData(targetAccumulator.RATIO.index)).getData()[1];
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
