/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.history;


/**
 * History that records a number of values, with new ones replacing the
 * earliest ones when the record is full.  The data returned by the
 * getHistory method will have the earliest first in the array, and the
 * most recent last.  If presented as a plot, the effect will be to 
 * scroll the data across the plot window. 
 * @author kofke, schultz, cribbin
 */
public class HistoryScrolling implements History, java.io.Serializable {
    
    public HistoryScrolling() {this(100);}
    public HistoryScrolling(int n) {
	    setHistoryLength(n);
	    reset();
    }
    
    /**
     * Sets the number of values kept in the history.
     */
    public void setHistoryLength(int n) {
    	int originalLength = history.length;
    	if (n==originalLength) return;
        
        // force current history into temp arrays
        getXValues();
        getHistory();
        
        xValues = new double[n];
        history = new double[n];
        if (n > originalLength) {
        	System.arraycopy(tempY,0,history,0,originalLength);
            System.arraycopy(tempX,0,xValues,0,originalLength);
        	for (int i = originalLength; i<n; i++) {
        		xValues[i] = Double.NaN;
                history[i] = Double.NaN;
        	}
        	cursor = originalLength;
        	full = false;
        }
        else {
        	System.arraycopy(tempX,originalLength-n,xValues,0,n);
            System.arraycopy(tempY,originalLength-n,history,0,n);
            if (cursor > n) {
                cursor = 0;
                full = true;
            }
        }
        tempX = new double[n];
        tempY = new double[n];
    }
    
    public int getHistoryLength() {
        return history.length;
    }
    
    public int getSampleCount() {
        return full ? history.length : cursor;
    }
	
    /**
     * Removes entire history, setting all values to NaN.
     */
	public void reset() {
	    for(int i=0; i<history.length; i++) {
	        history[i] = Double.NaN;
            xValues[i] = Double.NaN;
	    }
	    cursor = 0;
	    full = false;
	}
    
    public double[] getXValues() {
        int n=history.length;

        System.arraycopy(xValues,cursor,tempX,0,n-cursor);
        System.arraycopy(xValues,0,tempX,n-cursor,cursor);
        
        return tempX;
    }
	
    public boolean addValue(double x, double y) {
        xValues[cursor] = x;
        history[cursor] = y;
        cursor++;
        if(cursor == history.length) {
            cursor = 0;
            full = true;
        }
        return true;
    }

    /**
     * Returns an array with the most recent history at the end.
     */
    public double[] getHistory() {
		int n=history.length;

		System.arraycopy(history,cursor,tempY,0,n-cursor);
		System.arraycopy(history,0,tempY,n-cursor,cursor);
		
		return tempY;
    }
    
    public boolean willDiscardNextData() {
        return false;
    }

    private static final long serialVersionUID = 1L;
    private double[] history = new double[0];
    private int cursor;
	private double[] tempY = new double[0];
    private double[] tempX = new double[0];
    private double[] xValues = new double[0];
    private boolean full = false;

}
