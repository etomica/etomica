/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

import java.util.Arrays;

import static java.lang.Math.sqrt;

/**
 * Iterates arrays of int such that a[0] >= a[1] >= ... >= a[D-1] >= 0, where
 * D is the length of the array, which is set at construction.  Termination is
 * specified via the maxElement field.
 * For example, if D is 3, and maxElement is 2, iterator will return with successive calls to next():<br>
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
 *
 * The iterator may be set to disallow equal elements, so that iterates instead obey a[0] > a[1] > ... > a[D-1] >= 0.
 * Then, for example, if D = 3 and maxElement is 4, successive calls to next() return<br>
 * {2,1,0}<br>
 * {3,1,0}<br>
 * {3,2,0}<br>
 * {3,2,1}<br>
 * {4,1,0}<br>
 * {4,2,0}<br>
 * {4,2,1}<br>
 * {4,3,0}<br>
 * {4,3,1}<br>
 * {4,3,2}<br>
 *
 * @author David Kofke
 *
 */

/*
 * Note that IndexIteratorTriangular(D, maxElement, allowEqualElements=false) returns the same set of arrays as
 * CombinationIterator(n = maxElement+1, k = D) (but returned in a different order, and with the elements permuted)
 */ 
public class IndexIteratorTriangular implements IndexIterator, java.io.Serializable {

    /**
     * Constructs iterator that will return int arrays of length D.  Default value
     * of maxElement is D-1.
     */
    public IndexIteratorTriangular(int D) {
        if(D < 0) throw new IllegalArgumentException("Iterator must have non-negative dimension. Given value is "+D);
        this.D = D;
        index = new int[D];
        indexCopy = new int[D];
        setAllowEqualElements(true);
        setMaxElement(D-1);
        setMaxElementMin(0);
        reset();
    }

    /**
     * Indicates whether prespecified number of iterate have been returned
     * since last call to reset().
     */
    public boolean hasNext() {
        return hasNext;
    }
    
    /**
     * Returns the next iterate.  Status of hasNext is not relevant to behavior, and
     * will continue return new iterates with repeated calls.  The same array instance 
     * is returned with each call to this method; only its element values are modified in accordance
     * with the design of the iteration.
     */
    public int[] next() {
        increment(index, D-1);
        System.arraycopy(index, 0, indexCopy, 0, index.length);
        return indexCopy;
    }
    
    /**
     * Sets iterator to begin returning iterates.
     */
    public void reset() {
        if(allowEqualElements) {
            for(int i=0; i<D; i++) {
                index[i] = 0;
            }
            hasNext = maxElement >= maxElementMin;
        } else {
            for(int i=0; i<D; i++) {
                index[i] = D-1-i;
            }
            hasNext = (maxElement >= D-1) && (maxElement >= maxElementMin);
        }
        if(D == 1) {
            index[0] = Math.max(0, maxElementMin)-1;
        } else if(D > 1) {
            index[0] = Math.max(index[0], maxElementMin);
            index[D-1] = -1;
        } else {
            hasNext = false;
        }
    }
    
    private void increment(int[] idx, int d) {
        idx[d]++;
        if(allowEqualElements) {
            while(d > 0 && idx[d] > idx[d-1]) {
                idx[d] = 0;
                idx[--d]++;//decrement d, then increment idx
            }
            hasNext = (idx[D-1] < maxElement);
        } else {
            while(d > 0 && idx[d] >= idx[d-1]) {
                idx[--d]++;//decrement d, then increment idx
            }
            for(int i=d+1; i<D; i++) {
                idx[i] = D-1-i;
            }
            hasNext = (idx[0] < maxElement || idx[D-1] < maxElement-(D-1));
        }
    }

    /**
     * Returns the length of the int array that is returned as each iterate.
     */
    public int getD() {
        return D;
    }
 
    /**
     * Sets flag indicating whether iterates may have elements that are equal to one another.
     */
    public void setAllowEqualElements(boolean b) {
        allowEqualElements = b;
    }
    
    public boolean isAllowEqualElements() {
        return allowEqualElements;
    }

    public int getMaxElement() {
        return maxElement;
    }

    /**
     * Sets the largest allow value for any element to have.  If allowEqualElements,
     * the last iterate from the iterator will have all elements equal to this value;
     * if not allowEqualElements, the last iterate will have its zeroth element equal
     * to this value.
     */
    public void setMaxElement(int maxElement) {
        this.maxElement = maxElement;
    }

    public int getMaxElementMin() {
        return maxElementMin;
    }

    public void setMaxElementMin(int maxElementMin) {
        if(maxElementMin < 0) throw new IllegalArgumentException("Input value of maxElementMin must be nonnegative. Alas, this value was input: "+maxElementMin);
        this.maxElementMin = maxElementMin;
    }

    /**
     * Method to test and demonstrate class.
     */
    public static void main(String[] args) {
        IndexIteratorTriangular iterator = new IndexIteratorTriangular(3);
        iterator.setMaxElementMin(0);
        iterator.setMaxElement(4);
        //iterator.setAllowEqualElements(false);
        iterator.reset();
        System.out.println("Start");
        int count = 0;
        while(iterator.hasNext()) {
            int[] a = iterator.next();
            int sum = 0;
            for(int i=0; i<a.length; i++) {
                sum += a[i]*a[i];
            }
            System.out.println(++count+". "+ Arrays.toString(a) +" "+ sqrt(sum));
        }
    }

    private final int[] index;
    private final int[] indexCopy;
    private int D;
    private boolean allowEqualElements;
    private boolean hasNext;
    private int maxElement;
    private int maxElementMin;
    private static final long serialVersionUID = 1L;
}
