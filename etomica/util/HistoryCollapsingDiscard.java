package etomica.util;


/**
 * History that records a number of values.  When existing 
 * space is insufficient to hold new data, the existing data
 * is "collapsed" such that all existing data is stored in the first
 * half of existing storage by dropping every other data point.
 * After that point data will be taken half as often.
 * 
 * @author Andrew Schultz
 */
public class HistoryCollapsingDiscard extends HistoryCollapsing {

    public HistoryCollapsingDiscard() {this(100);}
    public HistoryCollapsingDiscard(int n) {
        this(n, 2);
    }

    public HistoryCollapsingDiscard(int nBins, int nCollapseBins) {
        super(nBins, nCollapseBins);
    }

    /**
     * adds data to the history.  If insufficient space exists
     * to store the data, existing data is collapsed by 1/2 and 
     * future data is taken half as often.
     */
    public void addValue(double x, double y) {
        if (++intervalCount == (interval+1)/2) {
            if (cursor == history.length) {
                collapseData();
            }
            xValues[cursor] = x;
            history[cursor] = y;
            cursor++;
        }
        if (intervalCount == interval) {
            intervalCount = 0;
        }
    }

    protected void collapseData() {
        for (int i=0; i<cursor/numCollapseBins; i++) {
            // j is the middle bin of the to-be-collapsed bins
            int j = i*numCollapseBins+(numCollapseBins-1)/2;
            history[i] = history[j];
            xValues[i] = xValues[j];
        }
        for(int i=cursor/numCollapseBins; i<cursor; i++) {
            history[i] = Double.NaN;
            xValues[i] = Double.NaN;
        }
        cursor /= numCollapseBins;
        interval *= numCollapseBins;
    }
}