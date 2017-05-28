/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.history;


/**
 * History that records a data from a simulation.  When existing space
 * is insufficient to hold new data, the array to hold the data is doubled
 * in length.
 * 
 * @author Andrew Schultz
 */
public class HistoryComplete implements History {
    
    public HistoryComplete() {this(100);}
    public HistoryComplete(int n) {
        setHistoryLength(n);
        reset();
        doCollapseOnReset = true;
    }
    
    /**
     * Sets the number of values kept in the history.
     */
    public void setHistoryLength(int n) {
    	int originalLength = history.length;
        if (n < cursor) {
            throw new IllegalArgumentException("decreasing history length would clobber data");
        }
    	if (n==originalLength) return;

        double[] temp = new double[n];
        System.arraycopy(xValues,0,temp,0,cursor);
        xValues = temp;
        for (int i = cursor; i<n; i++) {
            xValues[i] = Double.NaN;
        }

        temp = new double[n];
        System.arraycopy(history,0,temp,0,cursor);
        history = temp;
        for (int i = cursor; i<n; i++) {
            history[i] = Double.NaN;
        }
    }

    public int getHistoryLength() {
        return history.length;
    }
    
    public int getSampleCount() {
        return cursor;
    }
	
    /**
     * Removes entire history, setting all values to NaN.
     */
    public void reset() {
        int nValues = getHistoryLength();
        if (doCollapseOnReset && nValues > 100) {
            xValues = new double[100];
            history = new double[100];
        }
        
        for(int i=0; i<nValues; i++) {
            xValues[i] = Double.NaN;
            history[i] = Double.NaN;
        }
        cursor = 0;
    }
    
    public double[] getXValues() {
        return xValues;
    }
    
    public boolean addValue(double x, double y) {
        if (cursor == history.length) {
            int newLength = history.length*2;
            if (newLength == 0) {
                newLength = 1;
            }
            setHistoryLength(newLength);
        }
        xValues[cursor] = x;
        history[cursor] = y;
        cursor++;
        return true;
    }

    /**
     * Returns an the history
     */
    public double[] getHistory() {
        return history;
    }
    
    /**
     * Sets a flag that determines if the history is collapsed (to 100 data 
     * points) when reset is called.  If the current history length is less 
     * than 100, the history length is unchanged.
     */
    public void setCollapseOnReset(boolean flag) {
        doCollapseOnReset = flag;
    }
    
    /**
     * Returns true if the history is collapsed (to 100 data points) when reset
     * is called.  If the current history length is less than 100, the history 
     * length is unchanged.
     */
    public boolean isCollapseOnReset() {
        return doCollapseOnReset;
    }

    public boolean willDiscardNextData() {
        return false;
    }

    private double[] history = new double[0];
    private double[] xValues = new double[0];
    private int cursor;
    private boolean doCollapseOnReset;

}
