package etomica.util;

import etomica.data.DataSource;
import etomica.data.DataSourceUniform;

/**
 * History that records a number of values.  When existing 
 * space is insufficient to hold new data, the existing data
 * is "collapsed" such that all existing data is stored in the first
 * half of existing storage by dropping every other data point.
 * After that point data will be taken half as often.
 * 
 * @author Andrew Schultz
 */
public class HistoryCollapsing implements History {
    
    public HistoryCollapsing() {this(100);}
    public HistoryCollapsing(int n) {
        xSource = new DataSourceUniform();
        xSource.setTypeMin(DataSourceUniform.INCLUSIVE);
        xSource.setTypeMax(DataSourceUniform.INCLUSIVE);
        xSource.setXMin(0.0);
        setHistoryLength(n);
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
        if (n < 2) {
            throw new IllegalArgumentException("You have GOT to be kidding");
        }
        while (n < cursor) {
            collapseData();
        }
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
     * After calling reset, data will be recorded at every interval.
     */
    public void reset() {
        int nValues = getHistoryLength();
        for(int i=0; i<nValues; i++) {
            history[i] = Double.NaN;
        }
        cursor = 0;
        interval = 1;
        intervalCount = 0;
    }
    
    public DataSource getXSource() {
        return xSource;
    }
	
    /**
     * adds data to the history.  If insufficient space exists
     * to store the data, existing data is collapsed by 1/2 and 
     * future data is taken half as often.
     */
    public void addValue(double x) {
        if (++intervalCount != interval) {
            return;
        }
        if (cursor == history.length) {
            collapseData();
        }
        history[cursor] = x;
        cursor++;
    }

    protected void collapseData() {
        for (int i=0; i<cursor/2; i++) {
            history[i] = history[i*2];
        }
        cursor /= 2;
        interval *= 2;
    }
    
    /**
     * Returns the history
     */
    public double[] getHistory() {
        return history;
    }
    	
    /**
     * Factory that creates an instance of this class.
     */
    public static final History.Factory FACTORY = new History.Factory() {
        public History makeHistory() {return new HistoryCollapsing();}
        public History makeHistory(int n) {return new HistoryCollapsing(n);}
    };
	
    protected double[] history = new double[0];
    protected int cursor;
    private final DataSourceUniform xSource;
    protected int interval = 1;
    protected int intervalCount = 0;
}