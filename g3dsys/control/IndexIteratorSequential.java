package g3dsys.control;


/**
 * Simple sizeable index iterator that is in effect a set of nested loops.
 * index[D-1] iterates from 0 to size[D-1], and this is nested in iteration
 * over index[D-2], etc.  For example, if size = {3,2}, iteration would give
 * these values for index<br>
 * (0,0), (0,1), (1,0), (1,1), (2,0), (2,1).
 */

/*
 * History
 * Created on Jan 6, 2005 by kofke
 */
public class IndexIteratorSequential {

  
    /**
     * Constructs iterator that gives an index array of length D.
     * Default size is 1 for each dimension.  setSize and reset are
     * needed before beginning iteration. 
     */
    public IndexIteratorSequential(int D, int numLayers) {
        this.D = D;
        size = new int[D];
        for(int i=0; i<D; i++) size[i] = 1;
        setSize(numLayers*2+1);
        index = new int[D];
    }

    /* (non-Javadoc)
     * @see etomica.lattice.IndexIterator#reset()
     */
    public void reset() {
        for(int i=0; i<D; i++) {
            index[i] = 0;
        }
        index[D-1] = -1;
        count = 0;
    }

    /* (non-Javadoc)
     * @see etomica.lattice.IndexIterator#hasNext()
     */
    public boolean hasNext() {
        return count < maxCount;
    }

    /* (non-Javadoc)
     * @see etomica.lattice.IndexIterator#next()
     */
    public int[] next() {
        if(!hasNext()) return null;
        increment(index, D-1);
        count++;
        return index;
    }
    
    /**
     * Sets range of all indexes to the given value
     * @param size
     */
    private void setSize(int size) {
        maxCount = 1;
        for(int i=0; i<D; i++) {
            this.size[i] = size;
            maxCount *= size;
        }
    }
    
    public int getD() {
        return D;
    }
    
    private void increment(int[] idx, int d) {
        idx[d]++;
        while(idx[d] == size[d] && d > 0) {
            idx[d] = 0;
            idx[--d]++;//decrement d, then increment idx
        }
    }

    private final int[] size;
    private final int[] index;
    private final int D;
    //count is the number of iterates given since last reset
    //when this equals maxCount, iteration is done
    private int count;
    private int maxCount;
    private static final long serialVersionUID = 1L;
}
