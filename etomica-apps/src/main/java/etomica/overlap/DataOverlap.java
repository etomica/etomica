/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.overlap;

import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsingLog;
import etomica.data.DataSplitter;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.overlap.IntegratorOverlap.ReferenceFracSource;
import etomica.util.Debug;
import etomica.math.numerical.AkimaSpline;
import etomica.math.numerical.LinearFit;

/**
 * Utility class for overlap sampling.  This class retrieves data from
 * reference and target systems and calculates averages, ratios and
 * uncertainties as well as corresponding quantities for the natural log.
 * This class retrieves its information for a DataSourceOverlapAvg, which
 * actually understands how to retrieve data from the reference or target systems.
 * 
 * @author Andrew Schultz
 */
public class DataOverlap implements ReferenceFracSource {

    protected DataSourceOverlapAvg refAvg, targetAvg;
	protected AlphaSource alphaSource;
	protected AkimaSpline akima;
	protected double[] lnAlpha, err, lnAlphaDiff, lnRatio;
	protected double[] tmp = new double[2];
	protected boolean includeAntibiasErr = true;

	protected static final int REF = 0, TARGET = 1, RATIO = 2;

	public DataOverlap(DataSourceOverlapAvg refAvg, DataSourceOverlapAvg targetAvg, AlphaSource alphaSource) {
	    this.refAvg = refAvg;
	    this.targetAvg = targetAvg;
        
        lnAlpha = new double[0];
        err = new double[0];
        lnAlphaDiff = new double[0];
        lnRatio = new double[0];
        this.alphaSource = alphaSource;
	}
	
	public void setRefAvg(DataSourceOverlapAvg newRefAvg) {
	    refAvg = newRefAvg;
	}

	public void setTargetAvg(DataSourceOverlapAvg newTargetAvg) {
        targetAvg = newTargetAvg;
    }

	/**
	 * Returns an object capable of returning the different alpha values being
	 * used.
	 */
	public AlphaSource getAlphaSource() {
	    return alphaSource;
	}

	protected void populateArrays(int which, boolean doLog) {
        int numAlpha = alphaSource.getNumAlpha();
        if (lnAlpha.length != numAlpha) {
            lnAlpha = new double[numAlpha];
            err = new double[numAlpha];
            lnAlphaDiff = new double[numAlpha];
            lnRatio = new double[numAlpha];
        }
        double lnNRatio = 0;
        if (which == RATIO) {
            double nRatio = refAvg.getSampleCount()/targetAvg.getSampleCount();
            if (nRatio == 0 || Double.isInfinite(nRatio)) {
                for (int j=0; j<numAlpha; j++) {
                    lnAlpha[j] = Double.NaN;
                    err[j] = Double.NaN;
                    lnAlphaDiff[j] = Double.NaN;
                    lnRatio[j] = Double.NaN;
                }
                return;
            }
            lnNRatio = Math.log(nRatio);
        }

        for (int j=0; j<numAlpha; j++) {
            double avg, jerr;
            if (which == RATIO) {
                if (doLog) {
                    avg = Math.log(getAverage(j));
                    jerr = getLogError(j);
                    if (jerr < 0) {
                        throw new RuntimeException("oops");
                    }
                }
                else {
                    avg = getAverage(j);
                    jerr = getError(j);
                }
            }
            else {
                DataSourceOverlapAvg source = which == TARGET ? targetAvg : refAvg;
                if (doLog) {
                    avg = ((DataSourceOverlapLogAvg)source).getLogAverage(j);
                    jerr = ((DataSourceOverlapLogAvg)source).getLogError(j);
                }
                else {
                    avg = source.getAverage(j);
                    jerr = source.getError(j);
                }
            }

            if (doLog) {
                lnRatio[j] = avg;
                err[j] = jerr;
            }
            else {
                lnRatio[j] = Math.log(avg);
                err[j] = jerr/avg;
            }
            lnAlpha[j] = Math.log(alphaSource.getAlpha(j));
            lnAlphaDiff[j] = lnRatio[j] - lnAlpha[j] + lnNRatio;
            if (dumpit) System.out.println(lnAlpha[j]+" "+lnAlphaDiff[j]+" "+err[j]);
        }
	}
	public boolean dumpit;
	
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
        int numAlpha = alphaSource.getNumAlpha();
	    if (numAlpha == 1) {
	        tmp[0] = getAverage(0);
	        tmp[1] = getError(0);
	        return tmp;
	    }

	    populateArrays(RATIO, false);
	    
        double newAlpha = 0, newErr = 0;
        if (lnAlphaDiff[0] < 0) {
            // first new alpha is less than initial first alpha
            newAlpha = lnAlpha[0]; //Math.exp(lnAlphaDiff[0] + lnAlpha[0]);
            newErr = Double.NaN; //newAlpha*err[0];
        }
        else if (lnAlphaDiff[numAlpha-1] > 0) {
            newAlpha = lnAlpha[numAlpha-1]; //Math.exp(lnAlphaDiff[numAlpha-1] + lnAlpha[numAlpha-1]);
            newErr = Double.NaN; //newAlpha*err[numAlpha-1];
        }
        else if (numAlpha > 3) {
            if (akima == null) {
                akima = new AkimaSpline();
            }
            akima.setInputData(lnAlpha, lnAlphaDiff);
            double min = lnAlpha[0];
            double max = lnAlpha[numAlpha-1];
            double ymin = lnAlphaDiff[0];
            double ymax = lnAlphaDiff[numAlpha-1];
            double[] x = new double[1];
            x[0] = min+(max-min)/(ymax-ymin)*(-ymin);
            while (x[0] > min && x[0] < max) {
                double y = akima.doInterpolation(x)[0];
                if (y == 0) {
                    break;
                }
                if (y*min > 0) {
                    max = x[0];
                    ymax = y;
                }
                else {
                    min = x[0];
                    ymin = y;
                }
                x[0] = min+(max-min)/(ymax-ymin)*(-ymin);
            }
            akima.setInputData(lnAlpha, err);
            newErr = akima.doInterpolation(x)[0];
            newAlpha = x[0];
        }
        else {
            //linear interpolation (only 3 points)
            for (int i=0; i<numAlpha; i++) {
                if (lnAlphaDiff[i] > 0 && lnAlphaDiff[i+1] < 0) {
                    double ix = lnAlpha[i]+(lnAlpha[i+1]-lnAlpha[i])/(lnAlphaDiff[i+1]-lnAlphaDiff[i])*(-lnAlphaDiff[i]);
                    newAlpha = ix;
                    newErr = err[i] + (err[i+1]-err[i])*(ix-lnAlpha[i])/(lnAlpha[i+1]-lnAlpha[i]);
                }
            }
        }

        if (newAlpha < lnAlpha[0]) {
            newAlpha = lnAlpha[0];
        }
        else if (newAlpha > lnAlpha[numAlpha-1]) {
            newAlpha = lnAlpha[numAlpha-1];
        }
        double nRatio = refAvg.getSampleCount()/targetAvg.getSampleCount();
        tmp[0] = Math.exp(newAlpha) / nRatio;
        tmp[1] = tmp[0] / nRatio * newErr;
        return tmp;
	}

    /**
     * Returns the natural log of the average ratio and its uncertainty for the
     * given value of alpha
     */
    public double[] getLogAverageAndError(double iAlpha) {
        return getAverageAndError(RATIO, iAlpha, true);
    }
    
    /**
     * Returns the natural log of the average value from the reference or target
     * system (as specified) and its uncertainty for the given value of alpha
     */
    public double[] getLogAverageAndError(boolean isReference, double iAlpha) {
        return getAverageAndError(isReference ? REF : TARGET, iAlpha, true);
    }
    
    /**
     * Returns the average ratio and its uncertainty for the given value of
     * alpha
     */
    public double[] getAverageAndError(double iAlpha) {
        return getAverageAndError(RATIO, iAlpha, false);
    }
    
    /**
     * Returns the average from the reference or target system (as specified)
     * and its uncertainty for the given value of alpha
     */
    public double[] getAverageAndError(boolean isReference, double iAlpha) {
        return getAverageAndError(isReference ? REF : TARGET, iAlpha, false);
    }
    
	protected double[] getAverageAndError(int which, double iAlpha, boolean doLog) {

        populateArrays(which, doLog);

        int numAlpha = alphaSource.getNumAlpha();

        double newRatio = 0, newErr = 0;
        double lnNewAlpha = Math.log(iAlpha);
        if (numAlpha > 1) {
            if (lnNewAlpha < lnAlpha[0]) {
                lnNewAlpha = lnAlpha[0];
            }
            else if (lnNewAlpha > lnAlpha[numAlpha-1]) {
                lnNewAlpha = lnAlpha[numAlpha-1];
            }
        }
        if (numAlpha > 3) {
            if (akima == null) {
                akima = new AkimaSpline();
            }
            akima.setInputData(lnAlpha, lnRatio);
            newRatio = akima.doInterpolation(new double[]{lnNewAlpha})[0];
//            if (Double.isNaN(newRatio)) {
//                throw new RuntimeException("oops");
//            }
            akima.setInputData(lnAlpha, err);
            newErr = akima.doInterpolation(new double[]{lnNewAlpha})[0];
            if (newErr < 0) {
                // spline curve manufactured a minimum below 0
                // just value closest to our newAlpha
                System.err.println("error interpolation failed");
                for (int i=0; i<numAlpha; i++) {
                    if (lnAlpha[i] <= lnNewAlpha && lnAlpha[i+1] >= lnNewAlpha) {
                        if (lnNewAlpha - lnAlpha[i] < lnAlpha[i+1] - lnNewAlpha) {
                            newErr = err[i];
                        }
                        else {
                            newErr = err[i+1];
                        }
                        break;
                    }
                }
            }
        }
        else if (numAlpha > 1) {
            //linear interpolation (only 3 points)
            for (int i=0; i<numAlpha-1; i++) {
                if (lnAlpha[i] <= lnNewAlpha && lnAlpha[i+1] >= lnNewAlpha) {
                    double ix = lnRatio[i] + (lnRatio[i+1] - lnRatio[i]) / (lnAlpha[i+1] - lnAlpha[i]) * (lnNewAlpha - lnAlpha[i]);
                    newRatio = ix;
                    newErr = err[i] + (err[i+1] - err[i]) / (lnAlpha[i+1] - lnAlpha[i]) * (lnNewAlpha - lnAlpha[i]);
                }
            }
        }
        else {
            newRatio = lnRatio[0];
            newErr = err[0];
        }

        if (doLog) {
            tmp[0] = newRatio;
//            if (Double.isNaN(newRatio)) {
//                throw new RuntimeException("oops");
//            }
            if (newErr < 0) {
                throw new RuntimeException("oops");
            }
            tmp[1] = newErr;
        }
        else {
            tmp[0] = Math.exp(newRatio);
            tmp[1] = tmp[0] * newErr;
        }
        return tmp;
	}

	/**
	 * Returns the ideal fraction of time (or steps) spent in the reference
	 * system under the assumption that the current data while spending
	 * oldFrac amount of the calculation on the reference system.
	 */
    public double getIdealRefFraction(double oldFrac) {
        int numAlpha = alphaSource.getNumAlpha();
        double optimalAlpha;
        if (numAlpha == 1) {
            optimalAlpha = alphaSource.getAlpha(0);
        }
        else {
            optimalAlpha = getOverlapAverageAndError()[0];
        }
        if (Double.isNaN(optimalAlpha)) return 0.5;

        double refLnErr = getLogAverageAndError(true, optimalAlpha)[1];
        if (Double.isNaN(refLnErr)) {
            // if we don't have enough data to calc the error, assume it's 100%
            // if apparent error is > 100%, cap it there (>100% just means our estimate for ratio is bad)
            refLnErr = 1.0;
        }
        else if (includeAntibiasErr) {
            double refChi = getAverageAndError(true, optimalAlpha)[0];
            double refCount = refAvg.getSampleCount();
            double antibiasRefChi = (refChi * refCount + 1) / (refCount+1);
            double antibiasRefErr = 0.5*Math.log(antibiasRefChi/refChi);
            refLnErr += antibiasRefErr;
        }

        double targetLnErr = getLogAverageAndError(false, optimalAlpha)[1];
        if (Double.isNaN(targetLnErr)) {
            targetLnErr = 1.0;
        }
        else if (includeAntibiasErr) {
            double targetChi = getAverageAndError(false, optimalAlpha)[0];
            double targetCount = targetAvg.getSampleCount();
            double antibiasTargetChi = (targetChi * targetCount + 1/optimalAlpha) / (targetCount+1);
            double antibiasTargetErr = 0.5*Math.log(antibiasTargetChi/targetChi);
            targetLnErr += antibiasTargetErr;
        }
        double refFrac = 1.0 / (1 + targetLnErr/refLnErr * Math.sqrt((1-oldFrac)/(oldFrac)));
        return refFrac;
    }

    /**
     * Returns object capable of retrieving data from the reference system.
     */
    public DataSourceOverlapAvg getRefSource() {
        return refAvg;
    }

    /**
     * Returns object capable of retrieving data from the target system.
     */
    public DataSourceOverlapAvg getTargetSource() {
        return targetAvg;
    }

    /**
     * Returns the ratio of the reference to target overlap ratio
     * reference/target for the given value alpha.
     */
    public double getAverage(int iAlpha) {
        double t = targetAvg.getAverage(iAlpha);
        double r = refAvg.getAverage(iAlpha);
        return r/t;
    }

    /**
     * Returns the uncertainty in the log of the average corresponding to the
     * iAlpha value of alpha.
     */
    public double getLogError(int iAlpha) {
        if  (!(refAvg instanceof DataSourceOverlapAvgCollapsing)) {
            // just use standard propagation of error and pretend it works!
            return getError(iAlpha) / getAverage(iAlpha);
        }            
        double refLnErr = ((DataSourceOverlapAvgCollapsing)refAvg).getLogError(iAlpha);
        double targetLnErr = ((DataSourceOverlapAvgCollapsing)targetAvg).getLogError(iAlpha);
        return Math.sqrt(refLnErr*refLnErr + targetLnErr*targetLnErr);
    }
    
    /**
     * Returns the error in the ratio of the reference to target 
     * overlap-to-virial ratios (which reduces to target/reference) 
     * for the given value of the Bennet parameter.
     */
	public double getError(int iAlpha) {
		double refErr = refAvg.getError(iAlpha);
        double r = refAvg.getAverage(iAlpha);
        double refRelErr = refErr/r;
        double targetErr = targetAvg.getError(iAlpha);
        double t = targetAvg.getAverage(iAlpha);
        double targetRelErr = targetErr/t;
		return r/t*Math.sqrt(refRelErr*refRelErr+targetRelErr*targetRelErr);
	}

	/**
	 * Interface for a class capable of return the average and uncertainty for
	 * a particular value of alpha.
	 */
    public interface DataSourceOverlapAvg {
        public double getAverage(int iAlpha);
        public double getError(int iAlpha);
        public double getSampleCount();
    }
    
    /**
     * Interface for a class capable of return the average and uncertainty as
     * well as the log of the average and the uncertainty in the log, for
     * a particular value of alpha.
     */
    public interface DataSourceOverlapLogAvg extends DataSourceOverlapAvg {
        public double getLogError(int iAlpha);
        public double getLogAverage(int iAlpha);
    }
    
    /**
     * Implementation of DataSourceOverlapAvg that retrieves data from a single
     * AccumulatorAverage
     */
    public static class DataSourceOverlapAvgSimple implements DataSourceOverlapAvg {
        protected final AccumulatorAverage acc;
        public DataSourceOverlapAvgSimple(AccumulatorAverage acc) {
            this.acc = acc;
        }
        public double getAverage(int iAlpha) {
            return ((DataGroup)acc.getData()).getData(acc.AVERAGE.index).getValue(iAlpha);
        }
        public double getError(int iAlpha) {
            return ((DataGroup)acc.getData()).getData(acc.ERROR.index).getValue(iAlpha);
        }
        public double getSampleCount() {
            return acc.getBlockCount()*acc.getBlockSize();
        }
    }

    /**
     * Implementation of DataSourceOverlapLogAvg that retrieves data from a
     * splitter whose data is dumped into different
     * AccumulatorAverageCollapsingLog instances.  For each quantity, the
     * accumulator only provides results for a number of samples as a power of
     * 2.  This class attempts to estimate appropriate values for those 
     * quantities (primarily through extrapolation) to correspond to the total
     * number of samples collected.
     */
    public abstract static class DataSourceOverlapAvgCollapsing implements DataSourceOverlapLogAvg {
        protected final double[] lnStdev = new double[4];
        protected final double[] lnN = new double[4];
        protected double[] log = new double[0];
        protected double[] fullLnN = new double[0];
        protected final AkimaSpline akima = new AkimaSpline();

        protected abstract AccumulatorAverageCollapsingLog getDataSink(int i);

        protected abstract int getNumDataSinks();

        public double getAverage(int iAlpha) {
            IData data = getDataSink(iAlpha).getAverages();
            if (data.getLength() == 0) {
                return Double.NaN;
            }
            return data.getValue(0);
        }
        public double getError(int iAlpha) {
            AccumulatorAverageCollapsingLog acc = getDataSink(iAlpha);
            IData data = acc.getStdev();
            if (data.getLength() == 0) {
                return Double.NaN;
            }
            double count = acc.getCount();
            double count2 = 1L << (data.getLength()-1);
            return data.getValue(data.getLength()-1) * Math.sqrt(count2/count);
        }
        public double getLogError(int iAlpha) {
            return getExtrapolatedLogError(iAlpha, 0);
        }
        
        public double getExtrapolatedLogError(int iAlpha, int nDrop) {
            AccumulatorAverageCollapsingLog acc = getDataSink(iAlpha);
            IData data = acc.getStdevLog();
            int nData = data.getLength();
            if (nData == 0) {
                return Double.NaN;
            }
            int nFit = 4;
            if (nData <  nDrop + nFit || 1<<(nData-1) == acc.getCount()) {
                return data.getValue(nData-1);
            }
            boolean doDebug = Debug.ON && Math.random() > 0.98 && iAlpha == (getNumDataSinks()-1)/2;
            if (doDebug) {
                System.out.println("iAlpha stdev"+iAlpha+" "+nData);
                for (int i=0; i<nData; i++) {
                    System.out.println(i+" "+Math.log(1L<<i)+" "+Math.log(data.getValue(i)));
                }
            }
            double lnNFinal = Math.log(acc.getCount());
            int start = data.getLength() - nDrop - nFit;
            for (int i=0; i<nFit; i++) {
                lnN[i] = Math.log(1L << (start+i));
                lnStdev[i] = Math.log(data.getValue(start+i));
            }
            double[] fit = LinearFit.doFit(lnN, lnStdev);
            double m = fit[1];
            double b = fit[0];
            if (doDebug) {
                System.out.println("fit "+start+" "+m+" "+b);
            }
            if (m > 0) {
                if (doDebug) {
                    System.out.println("bogus m "+m+" "+b);
                }
                if (m < -0.5) {
                    m = -0.5;
                }
                else {
                    m = 0;
                }
                // refit with m fixed at new value
                double sum = 0;
                for (int i=0; i<nFit; i++) {
                    sum += lnStdev[i] - m*lnN[i];
                }
                b = sum / nFit;
            }
            double stdev = Math.exp(m*lnNFinal + b);
            if (doDebug) System.out.println(m+" "+b+" "+lnNFinal+" "+stdev+" "+data.getValue(data.getLength()-1));
            return stdev;
        }
        
        public double getLogAverage(int iAlpha) {
            AccumulatorAverageCollapsingLog acc = getDataSink(iAlpha);
            IData data = acc.getAverageLogs();
            int nData = data.getLength();
            if (nData == 0) {
                return Double.NaN;
            }
            if (nData <  11) {
                return data.getValue(nData-1);
            }
            boolean doDebug = Debug.ON && Math.random() > 0.99 && iAlpha == (getNumDataSinks()-1)/2;
            if (doDebug) {
                System.out.println("iAlpha log "+iAlpha+" "+nData);
                for (int i=0; i<nData; i++) {
                    System.out.println(i+" "+Math.log(1L<<i)+" "+data.getValue(i));
                }
            }

            // we fit our bias to A * exp(-B * l)
            //   where l = ln(N)
            
            if (fullLnN.length != data.getLength()) {
                fullLnN = new double[data.getLength()];
                log = new double[data.getLength()];
            }
            double ln2 = Math.log(2);
            for (int i=0; i<nData; i++) {
                fullLnN[i] = ln2*i;
                log[i] = data.getValue(i);
            }
            int nFit = 4;
            int nDrop = 5;
            int start = data.getLength() - nDrop - nFit;
            for (int i=0; i<nFit; i++) {
                lnN[i] = Math.log(1L << (start+i));
            }
            akima.setInputData(fullLnN, log);
            double[] dfdlnN = akima.doInterpolationDy(lnN);
            boolean isPositive = dfdlnN[0] > 0;
            if (doDebug) System.out.println("dfdlnN");
            for (int i=0; i<nFit; i++) {
                if (doDebug) {
                    System.out.println(fullLnN[start+i]+" "+dfdlnN[i]);
                }
                dfdlnN[i] = Math.log(Math.abs(dfdlnN[i]));
            }
            
            double[] fit = LinearFit.doFit(lnN, dfdlnN);
            double m = fit[1];
            double b = fit[0];
            if (doDebug) {
                System.out.println("fit "+start+" "+m+" "+b);
            }
            if (m < -1 || m > 0) {
                if (doDebug) {
                    System.out.println("bogus m "+m+" "+b);
                }
                if (m < -0.5) {
                    m = -0.5;
                }
                else {
                    // bias increasing with N, oops
                    m = 0.01;
                }
                // refit with m fixed at new value
                double sum = 0;
                for (int i=0; i<nFit; i++) {
                    sum += dfdlnN[i] - m*lnN[i];
                }
                b = sum / nFit;
            }
            double B = -m;
            double A = Math.exp(b)/B;
            if (!isPositive) A = -A;
            double avgLog = A*Math.exp(-B*lnN[nFit-1]) + log[start+nFit-1];
            if (doDebug) System.out.println(m+" "+b+" "+A+" "+avgLog+" "+data.getValue(data.getLength()-1));
//            if (Double.isNaN(avgLog)) {
//                throw new RuntimeException("oops");
//            }
            return avgLog;
        }

        public double getSampleCount() {
            return getDataSink(0).getCount();
        }
    }

    public static class DataSourceOverlapAvgCollapsingSplit extends DataSourceOverlapAvgCollapsing {
        protected final DataSplitter splitter;
        public DataSourceOverlapAvgCollapsingSplit(DataSplitter splitter) {
            this.splitter = splitter;
        }
        
        protected AccumulatorAverageCollapsingLog getDataSink(int i) {
            return (AccumulatorAverageCollapsingLog)splitter.getDataSink(i);
        }
        
        protected int getNumDataSinks() {
            return splitter.getNumDataSinks();
        }
    }
}
