/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.history;


/**
 * History that records a number of values.  When existing 
 * space is insufficient to hold new data, the existing data
 * is "collapsed" such that all existing data is stored in the first
 * half of existing storage by dropping every other data point.
 * After that point data will be taken half as often.
 * 
 * @author Andrew Schultz
 */
public abstract class HistoryCollapsing implements History {

    public HistoryCollapsing() {this(100);}
    public HistoryCollapsing(int n) {
        this(n, 2);
    }

    public HistoryCollapsing(int nBins, int nCollapseBins) {
        setNumCollapseBins(nCollapseBins);
        setHistoryLength(nBins);
        reset();
    }
    
    /**
     * Sets the number of values kept in the history.  If more data
     * is contained in the history than the set length, the existing 
     * data is collapsed until it fits within the new length.  History
     * length must be at least 2.
     */
    public void setHistoryLength(int n) {
        if (n==history.length) return;
        if (n < numCollapseBins) {
            throw new IllegalArgumentException("You have GOT to be kidding.  History length must be greater than "+numCollapseBins);
        }
        if (n % numCollapseBins != 0) {
            throw new IllegalArgumentException("History length ("+n+") must be an integer multiple of the # of collapse bins ("+numCollapseBins+")");
        }
        while (n < cursor) {
            collapseData();
        }

        double[] temp = new double[n];
        System.arraycopy(history,0,temp,0,cursor);
        history = temp;
        for (int i = cursor; i<n; i++) {
            history[i] = Double.NaN;
        }

        temp = new double[n];
        System.arraycopy(xValues,0,temp,0,cursor);
        xValues = temp;
        for (int i = cursor; i<n; i++) {
            xValues[i] = Double.NaN;
        }
    }

    public int getHistoryLength() {
        return history.length;
    }

    public int getSampleCount() {
        return cursor;
    }

    /**
     * Sets the number of "collapse" bins.  When the data is collapsed, one
     * value will be taken from each consecutive n bins, where n is the number
     * of bins set by this method.  The value will be taken from the middle bin
     * or (for an even number of bins), taken from the bin to the left of the
     * middle.
     */
    public void setNumCollapseBins(int newNumCollapseBins) {
        numCollapseBins = newNumCollapseBins;
        if (xValues.length % numCollapseBins != 0) {
            System.err.println("History length ("+xValues.length+") must be an integer multiple of the # of collapse bins ("+numCollapseBins+")");
        }
        reset();
    }
    
    public int getNumCollapseBins() {
        return numCollapseBins;
    }
    
    /**
     * Removes entire history, setting all values to NaN.
     * After calling reset, data will be recorded at every interval.
     */
    public void reset() {
        int nValues = getHistoryLength();
        for(int i=0; i<nValues; i++) {
            xValues[i] = Double.NaN;
            history[i] = Double.NaN;
        }
        cursor = 0;
        interval = 1;
        intervalCount = 0;
    }
    
    public double[] getXValues() {
        return xValues;
    }

    protected abstract void collapseData();

    /**
     * Returns the history
     */
    public double[] getHistory() {
        return history;
    }

    protected double[] history = new double[0];
    protected int cursor;
    protected double[] xValues = new double[0];
    protected int interval;
    protected int intervalCount;
    protected int numCollapseBins;
}
