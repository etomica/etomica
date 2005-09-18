package etomica.util;

import etomica.data.DataSource;
import etomica.data.DataSourceUniform;

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
        xSource = new DataSourceUniform();
        xSource.setTypeMin(DataSourceUniform.INCLUSIVE);
        xSource.setTypeMax(DataSourceUniform.INCLUSIVE);
        xSource.setXMin(0.0);
        setHistoryLength(n);
        reset();
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
    	xSource.setNValues(n);
        xSource.setXMax(n);
        double[] temp = new double[n];
        System.arraycopy(history,0,temp,0,cursor);
        history = temp;
        for (int i = cursor; i<n; i++) {
            history[i] = Double.NaN;
        }
    }

    public int getHistoryLength() {return history.length;}
	
    /**
     * Removes entire history, setting all values to NaN.
     */
    public void reset() {
        int nValues = getHistoryLength();
        for(int i=0; i<nValues; i++) {
            history[i] = Double.NaN;
        }
        cursor = 0;
    }
    
    public DataSource getXSource() {
        return xSource;
    }
    
    public void addValue(double x) {
        if (cursor == history.length) {
            int newLength = history.length*2;
            if (newLength == 0) {
                newLength = 1;
            }
            setHistoryLength(newLength);
        }
        history[cursor] = x;
        cursor++;
    }

    /**
     * Returns an the history
     */
    public double[] getHistory() {
        return history;
    }
    	
    /**
     * Factory that creates an instance of this class.
     */
    public static final History.Factory FACTORY = new History.Factory() {
        public History makeHistory() {return new HistoryComplete();}
        public History makeHistory(int n) {return new HistoryComplete(n);}
    };
    
    private double[] history = new double[0];
    private int cursor;
    private final DataSourceUniform xSource;

}