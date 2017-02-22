/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

import java.util.Arrays;


/**
 * Iterator that returns indices in increasing order, starting at 0 and ending
 * including all indices where |i| <= max, where i are the components of the
 * index.
 *
 * This class has the same effect as
 * IndexIteratorReflector + IndexIteratorTriangularPermutations, but is ~15x
 * faster.
 * 
 * @author Andrew Schultz
 */
public class IndexIteratorIncreasing {

    public IndexIteratorIncreasing(int dim) {
        this.dim = dim;
        next = new int[dim];
        prev = new int[dim];
    }
    
    /**
     * Resets the iterator.  The first iterate is always 0,0,0...
     */
    public void reset() {
        for (int i=0; i<dim; i++) {
            next[i] = 0;
            prev[i] = 0;
        }
        m = -1;
        firstM = 0;
    }
    
    /**
     * Returns the next index.  Returns null if there are no more indices.
     */
    public int[] next() {
        if (m==-1) {
            m = 0;
            return next;
        }
        for (int i = dim-1; i>=0; i--) {
            if (i==firstM) continue;
            if (next[i] < m && (i>firstM || next[i]<m-1)) {
                next[i]++;
                return next;
            }
            if (i>firstM) {
                next[i] = -m;
            }
            else {
                next[i] = -m+1;
            }
        }
        if (m > 0 && next[firstM] == -m) {
            next[firstM] = m;
            return next;
        }
        if (firstM < dim-1 && m > 0) {
            next[firstM] = -m+1;
            firstM++;
            return next;
        }
        m++;
        firstM = 0;
        if (m > max) {
            return null;
        }
        for (int i=0; i<dim; i++) {
            next[i] = -m;
        }
        return next;
    }

    /**
     * Sets a new maximum for iteration.  The iterator will return all sets of
     * indices such that |i|<m (where i is an index in the set).
     */
    public void setMax(int newMax) {
        if (newMax < 0) {
            throw new RuntimeException("must be positive");
        }
        max = newMax;
    }
    
    protected final int dim;
    protected final int[] next;
    protected final int[] prev;
    protected int max;
    protected int m, firstM;
    
    public static void main(String[] args) {
        long t1 = System.currentTimeMillis();
        IndexIteratorIncreasing iter = new IndexIteratorIncreasing(3);
        iter.setMax(Integer.parseInt(args[0]));
        iter.reset();
        int[] next;
        while ((next = iter.next()) != null) {
            System.out.println(Arrays.toString(next));
        }
        System.out.println(System.currentTimeMillis() - t1);
    }
        
}
