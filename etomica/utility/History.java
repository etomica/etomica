package etomica.utility;

import etomica.*;
import etomica.units.Dimension;

public class History implements DataSource {
    
    private double[] values;
    private int cursor;
    private String label;
    private Dimension dimension;
    private String name;
    
    public History() {this(100);}
    public History(int n) {
	    setNValues(n);
	    reset();
	    setDimension(Dimension.NULL);
    }
    
	public void setName(String s) {name = s;}
	public String toString() {return name;}
	
    public void setNValues(int n) {
        values = new double[n];
        cursor = 0;
    }
    public int getNValues() {return (values!=null) ? values.length : 0;}
	
	public void reset() {
	    int nValues = getNValues();
	    for(int i=0; i<nValues; i++) {
	        values[i] = Double.NaN;
	    }
	    cursor = 0;
	}
	
    public void addValue(double x) {
        values[cursor] = x;
        cursor++;
        if(cursor == values.length) cursor = 0;
    }
    
/*	 stuff copied from Accumulator's implementation
       public void setHistoryWindow() {
	        double[] oldHistory = history;
	        history = new double[historyWindow];
	        clearHistory(); //set to NaN
	        if(oldHistory != null) {
    	        int n = Math.min(oldHistory.length, historyWindow); //setHistoryWindow in MeterAbstract will not let this be less than 1
	            for(int i=0; i<n; i++) {
	                history[i] = oldHistory[i];
	            }
	            historyCursor = n;
	        }
	        else historyCursor = 0;
	    }//end of setHistoryWindow
	    
	            
        public void addToHistory(double value) {
            if(historyCursor < historyWindow) {
                history[historyCursor++] = value;
            }
            else {
                MeterAbstract.this.setHistoryWindow(2*historyWindow);
                addToHistory(value);
            }
        }
        
	/**
	 * Accessor method for the period of the history, i.e., the number of values
	 * that are kept recorded.
	 
	public int getHistoryWindow() {return historyWindow;}

	/**
	 * Accessor method for the period of the history, i.e., the number of values
	 * that are kept recorded.  Default is 100, minimum is 1.
	 
	public void setHistoryWindow(int window) {
	    if(window < 1) window = 1; 
	    historyWindow = window;
	    if(!isHistorying()) return;
	    Accumulator[] accumulators = allAccumulators();
	    for(int i=0; i<accumulators.length; i++) {
	        accumulators[i].setHistoryWindow();
	    }
	}
        
*/
        public double[] values(DataSource.ValueType dummy) {
            return values;
        }
        
        public etomica.DataSource.ValueType[] dataChoices() {return null;}
        
        public String getLabel() {return label;}
        public void setLabel(String text) {label = text;}
        
        public Dimension getDimension() {return dimension;}
        public void setDimension(Dimension dim) {dimension = dim;}
        
}
