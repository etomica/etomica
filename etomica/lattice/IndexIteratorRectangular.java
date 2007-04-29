package etomica.lattice;

import etomica.util.Arrays;

/**
 * Iterates arrays of int, such that 0 <= a[0] <= size[0], 0 <= a[1] <= size[1], etc., 
 * thereby spanning a rectangular region. For example, if size = {3,2}, iteration would give
 * these values and then expire:<br>
 * {0,0}<br>
 * {0,1}<br>
 * {1,0}<br>
 * {1,1}<br>
 * {2,0}<br>
 * {2,1}
 */

public class IndexIteratorRectangular implements IndexIterator, IndexIteratorSizable, java.io.Serializable {

    /**
     * Constructs iterator that returns int arrays of length D.
     * Default size is 1 for each dimension, so that only the iterate {0,0,..,0}
     * is returned; setSize can be used to modify the range.
     * Call to reset is required before beginning iteration, in any case. 
     */
    public IndexIteratorRectangular(int D) {
        this.D = D;
        size = new int[D];
        for(int i=0; i<D; i++) size[i] = 1;
        setSize(size);
        index = new int[D];
    }

    /**
     * Puts iterator in a condition to begin iteration.  First iterate is {0,...,0}
     */
    public void reset() {
        for(int i=0; i<D; i++) {
            index[i] = 0;
        }
        index[D-1] = -1;
        count = 0;
    }

    /**
     * Indicates whether all iterates have been returned since last call to reset().
     */
    public boolean hasNext() {
        return count < maxCount;
    }

    /**
     * Returns the next iterate, or null if hasNext is false. The same array instance 
     * is returned with each call to this method; only its element values are modified in accordance
     * with the design of the iteration.
     */
    public int[] next() {
        if(!hasNext()) return null;
        increment(index, D-1);
        count++;
        return index;
    }
    
    /**
     * Sets range of all indexes to the given value
     * Indices returned vary from 0 to size-1
     */
    public void setSize(int size) {
        maxCount = 1;
        for(int i=0; i<D; i++) {
            this.size[i] = size;
            maxCount *= size;
        }
    }
    
    /**
     * Sets the individual ranges for the elements of the returned
     * array.  Each element index[k] iterates over values from 0 to size[k]-1. 
     */
    public void setSize(int[] size) {
        if(size.length != D) throw new IllegalArgumentException("Length of array inconsistent with dimension of iterator");
        maxCount = 1;
        for(int i=0; i<D; i++) {
            this.size[i] = size[i];
            maxCount *= size[i];
        }
    }

    /**
     * Returns the length of the integer array that is returned as an iterate.
     */
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

    /**
     * main method to test and demonstrate class.
     */
    public static void main(String[] args) {
        IndexIteratorRectangular iterator = new IndexIteratorRectangular(3);
        iterator.reset();
        System.out.println("Start 1");
        while(iterator.hasNext()) {
            System.out.println(Arrays.toString(iterator.next()));
        }
        iterator.setSize(new int[] {5,2,4});
        iterator.reset();
        System.out.println("Start 2");
        while(iterator.hasNext()) {
            System.out.println(Arrays.toString(iterator.next()));
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
