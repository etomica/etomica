package etomica.util;


/**
 * History that records a number of values.  When existing 
 * space is insufficient to hold new data, the existing data
 * is "collapsed" such that all existing data is stored in the first
 * half of existing storage by dropping every other data point.
 * After that point data will be taken half as much.
 * 
 * @author Andrew Schultz
 */
public class HistoryCollapsingAverage extends HistoryCollapsing {
    
    public HistoryCollapsingAverage() {this(100);}
    public HistoryCollapsingAverage(int n) {
        super(n);
    }
    
   /**
     * adds data to the history.  If insufficient space exists
     * to store the data, existing data is collapsed by 1/2 and 
     * future data is taken half as much.
     */
    public void addValue(double x) {
        if (cursor == history.length) {
            collapseData();
            history[cursor] = 0.0;
        }
        history[cursor] += x;
        if (++intervalCount == interval) {
            history[cursor] /= interval;
            cursor++;
            if (cursor < history.length) {
                history[cursor] = 0.0;
            }
        }
    }

    protected void collapseData() {
        for (int i=0; i<cursor/2; i++) {
            history[i] = (history[i*2] + history[i*2+1]);
        }
        for (int i=cursor/2; i<history.length; i++) {
            history[i] = Double.NaN;
        }
        cursor /= 2;
        interval *= 2;
    }

    /**
     * Returns the history
     */
    public double[] getHistory() {
        if (intervalCount == 0) {
            return history;
        }
        if (temp == null || temp.length != history.length) {
            temp = new double[history.length];
        }
        System.arraycopy(history,0,temp,0,history.length);
        history[cursor] /= intervalCount;
        return temp;
    }
    
    protected double[] temp;
}