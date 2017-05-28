/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.histogram;

import etomica.math.DoubleRange;

/**
 * Histogram class capable of keeping track of how many times each discrete
 * value was passed in, similar to what one would expect from a histogram with
 * very small bin width.  Two incoming values are judged to be the same if they
 * are within some tolerance, which is specified at construction.  The returned
 * histogram is then the (raw) probability with which each value was seen.
 *   
 * @author Andrew Schultz
 */
public class HistogramDiscrete implements Histogram {

    public HistogramDiscrete(double tolerance) {
        this.tolerance = tolerance;
        histogram = new double[0];
        xValues = new double[0];
        counts = new long[0];
    }
    
    public void addValue(double x) {
        int bin;
        if (counts.length == 0) {
            counts = new long[1];
            xValues = new double[1];
            bin = 0;
            xValues[bin] = x;
        }
        else {
            int maxBin = counts.length-1;
            int minBin = 0;
            bin = (minBin+maxBin)/2;
            double diff = x - xValues[bin];
            if (Double.isNaN(diff) && !Double.isNaN(x)) {
                diff = 0;
            }
            while (Math.abs(diff) > tolerance) {
                if (diff>0) {
                    if (bin == maxBin) {
                        long[] newCounts = new long[counts.length+1];
                        System.arraycopy(counts, 0, newCounts, 0, maxBin+1);
                        if (maxBin < counts.length-1) {
                            System.arraycopy(counts, maxBin+1, newCounts, maxBin+2, counts.length-1-maxBin);
                        }
                        counts = newCounts;
                        double[] newXValues = new double[xValues.length+1];
                        System.arraycopy(xValues, 0, newXValues, 0, maxBin+1);
                        if (maxBin < xValues.length-1) {
                            System.arraycopy(xValues, maxBin+1, newXValues, maxBin+2, xValues.length-1-maxBin);
                        }
                        xValues = newXValues;
                        bin++;
                        xValues[bin] = x;
                        break;
                    }
                    minBin = bin+1;
                    bin = (minBin+maxBin)/2;
                    diff = x - xValues[bin];
                }
                else {
                    if (bin == minBin) {
                        long[] newCounts = new long[counts.length+1];
                        if (minBin>0) {
                            System.arraycopy(counts, 0, newCounts, 0, minBin);
                        }
                        System.arraycopy(counts, minBin, newCounts, minBin+1, counts.length-minBin);
                        counts = newCounts;
                        double[] newXValues = new double[xValues.length+1];
                        if (minBin>0) {
                            System.arraycopy(xValues, 0, newXValues, 0, minBin);
                        }
                        System.arraycopy(xValues, minBin, newXValues, minBin+1, xValues.length-minBin);
                        xValues = newXValues;
                        xValues[bin] = x;
                        break;
                    }
                    maxBin = bin-1;
                    bin = (minBin+maxBin)/2;
                    diff = x - xValues[bin];
                }
            }
        }
        counts[bin]++;
        
        totalCount++;
    }

    public long getCount() {
        return totalCount;
    }

    public double[] getHistogram() {
        if (histogram.length != counts.length) {
            histogram = new double[counts.length];
        }
        for (int i=0; i<counts.length; i++) {
            histogram[i] = ((double)counts[i])/totalCount;
        }
        return histogram;
    }

    public int getNBins() {
        return counts.length;
    }

    public DoubleRange getXRange() {
        if (xValues.length == 0) return null;
        return new DoubleRange(xValues[0], xValues[xValues.length-1]);
    }

    public void reset() {
        xValues = new double[0];
        counts = new long[0];
        histogram = new double[0];
        totalCount = 0;
    }

    public void setNBins(int n) {
        // hi-larious
    }

    public void setXRange(DoubleRange range) {
        // very funny
    }

    public double[] xValues() {
        return xValues;
    }

    protected long totalCount;
    protected double tolerance;
    protected double[] xValues;
    protected long[] counts;
    protected double[] histogram;
}
