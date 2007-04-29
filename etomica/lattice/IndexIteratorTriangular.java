package etomica.lattice;

import etomica.util.Arrays;


/**
 * Returns arrays of int such that int[0] >= int[1] >= ... >= int[D-1], where
 * D is the length of the array, which is set at construction.  For example, 
 * if D is 3, will return, with successive calls to next():
 * {0,0,0}<br>
 * {1,0,0}<br>
 * {1,1,0}<br>
 * {1,1,1}<br>
 * {2,0,0}<br>
 * {2,1,0}<br>
 * {2,1,1}<br>
 * {2,2,0}<br>
 * {2,2,1}<br>
 * {2,2,2}<br>
 * {3,0,0}<br>
 * etc.
 *
 * @author David Kofke
 *
 */

public class IndexIteratorTriangular implements IndexIterator, java.io.Serializable {

    /**
     * 
     */
    public IndexIteratorTriangular(int D) {
        this.D = D;
        setNumberOfIterates(Integer.MAX_VALUE);
        reset();
    }

    /**
     * Indicates whether prespecified number of iterate have been returned
     * since last call to reset().
     */
    public boolean hasNext() {
        return iterateCount < numberOfIterates;
    }
    
    /**
     * Returns the next iterate.  Status of hasNext is not relevant to behavior, and
     * will continue return new iterates with repeated calls.
     */
    public int[] next() {
        increment(index, D-1);
        iterateCount++;
        return index;
    }
    
    /**
     * Sets iterator to begin returning iterates, beginning with {0,...,0}, and sets
     * iterate count to zero.
     */
    public void reset() {
        index = new int[D];
        iterateCount = 0;
        index[D-1] = -1;
    }
    
    private void increment(int[] idx, int d) {
        idx[d]++;
        while(d > 0 && idx[d] > idx[d-1]) {
            idx[d] = 0;
            idx[--d]++;//decrement d, then increment idx
        }
    }

    /**
     * Returns the length of the int array that is returned as each iterate.
     */
    public int getD() {
        return D;
    }
 
    public int getNumberOfIterates() {
        return numberOfIterates;
    }

    /**
     * Sets the value of the number of iterates that can be returned before hasNext
     * begins to return false.
     */
    public void setNumberOfIterates(int numberOfIterates) {
        this.numberOfIterates = numberOfIterates;
    }

    /**
     * Method to test and demonstrate class.
     */
    public static void main(String[] args) {
        IndexIteratorTriangular iterator = new IndexIteratorTriangular(3);
        iterator.setNumberOfIterates(30);
        iterator.reset();
        System.out.println("Start 2");
        while(iterator.hasNext()) {
            int[] a = iterator.next();
            int sum = 0;
            for(int i=0; i<a.length; i++) {
                sum += a[i]*a[i];
            }
            System.out.println(Arrays.toString(a)+" "+Math.sqrt(sum));
        }
    }

    private int[] index;
    private int D;
    private int numberOfIterates;
    private int iterateCount;
    private static final long serialVersionUID = 1L;
}
