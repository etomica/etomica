/* History
 * 08/09/04 (DAK, AS, NC) new from History (when it was defined as a class)
 */
package etomica.utility;

/**
 * History that records a number of values, with new ones replacing the
 * earliest ones when the record is full.  The data returned by the
 * getHistory method will have the earliest first in the array, and the
 * most recent last.  If presented as a plot, the effect will be to 
 * scroll the data across the plot window. 
 * @author kofke, schultz, cribbin
 */
public class HistoryScrolling implements History {
    
    public HistoryScrolling() {this(100);}
    public HistoryScrolling(int n) {
	    setHistoryLength(n);
	    reset();
    }
    
    /**
     * Sets the number of values kept in the history.
     */
    public void setHistoryLength(int n) {
    	int originalLength = getHistoryLength();
    	if (n==originalLength) return;
    	xTemp = new double[n];
    	for (int i=0; i<n; i++) {
    		xTemp[i] = Double.NaN;
    	}
    	// calling getHistory() fills the temp array
    	getHistory();
        history = new double[n];
        if (n > originalLength) {
        	System.arraycopy(temp,0,history,0,originalLength);
        	for (int i = originalLength; i<n; i++) {
        		history[i] = Double.NaN;
        	}
        	cursor = originalLength;
        }
        else {
        	System.arraycopy(temp,originalLength-n,history,0,n);
        	cursor = 0;
        }
        temp = new double[n];
    }
    public int getHistoryLength() {return history.length;}
	
    /**
     * Removes entire history, setting all values to NaN.
     */
	public void reset() {
	    int nValues = getHistoryLength();
	    for(int i=0; i<nValues; i++) {
	        history[i] = Double.NaN;
	        xTemp[i] = Double.NaN;
	    }
	    cursor = 0;
	}
	
    public void addValue(double x) {
        history[cursor] = x;
        cursor++;
        count++;
        if(cursor == history.length) cursor = 0;
    }

    /**
     * Returns an array with the most recent history at the end.
     */
    public double[] getHistory() {
		int n=history.length;

		System.arraycopy(history,cursor,temp,0,n-cursor);
		System.arraycopy(history,0,temp,n-cursor,cursor);
		
		return temp;
    }
    
	public double[] xValues() {
		int nValues = getHistoryLength();
		xTemp[nValues-1] = count;
		for (int i=nValues-2; i>-1 && xTemp[i+1]>1; i--) {
			xTemp[i] = xTemp[i+1] - 1;
		}
		return xTemp;
	}
	
	/**
	 * Factory that creates an instance of this class.
	 */
    public static final History.Factory FACTORY = 
    	new History.Factory() {
    		public History makeHistory() {return new HistoryScrolling();}
    		public History makeHistory(int n) {return new HistoryScrolling(n);}
    	};
	
    private double[] history = new double[0];
    private int cursor;
	private double[] temp = new double[0];
	private double[] xTemp;
	private int count=0;

}