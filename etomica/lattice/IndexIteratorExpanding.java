package etomica.lattice;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jan 5, 2005 by kofke
 */
public class IndexIteratorExpanding implements IndexIterator, java.io.Serializable {

    /**
     * 
     */
    public IndexIteratorExpanding(int D) {
        this.D = D;
        ratios = new double[D];
        for(int i=0; i<D; i++) ratios[i] = 1.0;
        reset();
    }

    /* (non-Javadoc)
     * @see etomica.lattice.IndexIterator#hasNext()
     */
    public boolean hasNext() {
        return true;
    }
    /* (non-Javadoc)
     * @see etomica.lattice.IndexIterator#next()
     */
    public int[] next() {
        increment(index, D-1);
        return index;
    }
    /* (non-Javadoc)
     * @see etomica.lattice.IndexIterator#reset()
     */
    public void reset() {
        index = new int[D];
    }
    
    private void increment(int[] idx, int d) {
        idx[d]++;
        while(idx[d] == size[d] && d > 0) {//replaces recursive call
            idx[d] = 0;
            idx[--d]++;//decrement d, then increment idx
        }
    }

    private int[] index;
    private double[] ratios;
    private int[] size;
    private int D;
}
