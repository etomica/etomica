package etomica.utility;

public class HistoryScrolling implements History {
    
    private double[] values;
    private int cursor;
	private double[] temp;
	private double[] xTemp;
	private int count=0;
	
    public HistoryScrolling() {this(100);}
    public HistoryScrolling(int n) {
	    setNValues(n);
	    reset();
    }
    
    /**
     * sets the number of values kept in the history.
     */
    public void setNValues(int n) {
    	int originalLength = getNValues();
    	if (n==originalLength) return;
    	xTemp = new double[n];
    	for (int i=0; i<n; i++) {
    		xTemp[i] = Double.NaN;
    	}
    	// calling getHistory() fills the temp array
    	getHistory();
        values = new double[n];
        if (n > originalLength) {
        	System.arraycopy(temp,0,values,0,originalLength);
        	for (int i = originalLength; i<n; i++) {
        		values[i] = Double.NaN;
        	}
        	cursor = originalLength;
        }
        else {
        	System.arraycopy(temp,originalLength-n,values,0,n);
        	cursor = 0;
        }
        temp = new double[n];
    }
    public int getNValues() {return values.length;}
	
	public void reset() {
	    int nValues = getNValues();
	    for(int i=0; i<nValues; i++) {
	        values[i] = Double.NaN;
	        xTemp[i] = Double.NaN;
	    }
	    cursor = 0;
	}
	
    public void addValue(double x) {
        values[cursor] = x;
        cursor++;
        count++;
        if(cursor == values.length) cursor = 0;
    }

    /**
     * returns an array with the most recent values at the end of
     * the array
     */
    public double[] getHistory() {
		int n=values.length;

		System.arraycopy(values,cursor,temp,0,n-cursor);
		System.arraycopy(values,0,temp,n-cursor,cursor);
		
		return temp;
    }
    
	public double[] xValues() {
		int nValues = getNValues();
		xTemp[nValues-1] = count;
		for (int i=nValues-2; i>-1 && xTemp[i+1]>1; i--) {
			xTemp[i] = xTemp[i+1] - 1;
		}
		return xTemp;
	}
	
}