/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.discrete;

import etomica.lattice.IndexIterator;
import etomica.lattice.IndexIteratorTriangular;
import etomica.math.SpecialFunctions;
import java.util.Arrays;

/**
 * Iterator that returns unique permutations of an array of integers, with
 * degeneracy such that some of the integers may appear more than once in the array. 
 * Size of array is declared at construction. The setDegeneracy method then describes
 * the number of degenerate elements that should be present in each iterate. For
 * example, if constructed with n = 5 and the degeneracy is {2,3}, then
 * successive calls to next() return the arrays:<br>
 * {1,1,1,0,0}<br>
 * {1,1,0,1,0}<br>
 * {1,0,1,1,0}<br>
 * {0,1,1,1,0}<br>
 * {1,1,0,0,1}<br>
 * {1,0,1,0,1}<br>
 * {0,1,1,0,1}<br>
 * {1,0,0,1,1}<br>
 * {0,1,0,1,1}<br>
 * {0,0,1,1,1}<br>
 * 
 * The total number of iterates returned equals n!/product{degeneracy[i]!}
 * 
 * @author David Kofke
 */

/*
 * Description of algorithm: To generate int[] having n elements, this instance
 * will use a subiterator instance to generate permutations of n-d elements,
 * where d is the degeneracy of the element associated with this instance. For
 * each such permutation, this instance will insert its element at d positions
 * distributed through the 0, 1, 2,...,n-1 in the sequence returned by the
 * subiterator.
 */
public class PermutationDegenerateIterator implements IndexIterator, java.io.Serializable {

    /**
     * Constructs instance whose iterates are int[] of length n. Default
     * degeneracy is defined to have no repeated elements. User must call
     * reset() before starting iteration.
     */
    public PermutationDegenerateIterator(int n) {
        super();
        if (n < 0)
            throw new IllegalArgumentException(
                    "PermutationIterator constructor requires non-negative argument");
        this.n = n;
        nextIndex = new int[n];
        zero = new int[n];
        int[] degeneracy = new int[n];
        for (int i = 0; i < n; i++) {
            degeneracy[i] = 1;
        }
        setDegeneracy(degeneracy);
        hasNext = false;
    }

    // constructor for use by NullSinglet inner class
    private PermutationDegenerateIterator() {
        subPermutation = null;
        nextIndex = null;
        zero = null;
        n = 0;
    }
    
    /**
     * Invokes setDegeneracy(degeneracy, degeneracy.length), thus using all elements of given array.
     */
    public void setDegeneracy(int[] degeneracy) {
        setDegeneracy(degeneracy, degeneracy.length);
    }

    /**
     * Sets parameter that describes number of duplicate elements in the
     * iterates. For example, if degeneracy is {1,3,2}, then each iterate will
     * have 1 element that is 0 (because degeneracy[0]=1), 3 elements that are 1
     * (because degeneracy[1] = 3), and 2 elements that are 2 (because
     * degeneracy[2] = 2). Thus the array {1,0,1,2,2,1} would for this example
     * be one of the iterates eventually returned by next(). Must be followed by
     * a call to reset() before this method has any effect.
     *  
     * Degeneracy is taken as the first "length" elements of the given array.
     * Any elements degeneracy[i] for {@code i>=length} are ignored.
     * 
     * @throws IllegalArgumentException
     *             if the sum of the elements in degeneracy is not equal to the
     *             value of n set a construction.
     */
    public void setDegeneracy(int[] degeneracy, int length) {
        this.degeneracy = degeneracy.clone();
        int sum = 0;
        for (int i = 0; i < length; i++) {
            if (degeneracy[i] < 0)
                throw new IllegalArgumentException(
                        "All degeneracy values must be nonnegative. You gave these: "
                                + Arrays.toString(degeneracy));
            sum += degeneracy[i];
        }
        if (sum != n)
            throw new IllegalArgumentException(
                    "Degeneracy inconsistent with length of iterates.  Length of iterates: "
                            + n + "; degeneracy: "
                            + Arrays.toString(degeneracy) + ", using first "+length+" elements, which sums to "
                            + sum);
        element = length - 1;
        if (length == 0) {
            subPermutation = new NullSinglet();
        } else {
            int thisDegeneracy = degeneracy[length - 1];
            int[] subDegeneracy = new int[length - 1];
            System.arraycopy(degeneracy, 0, subDegeneracy, 0, length - 1);
            indexIterator = new IndexIteratorTriangular(thisDegeneracy);
            indexIterator.setAllowEqualElements(false);
            indexIterator.setMaxElement(n - 1);
            subPermutation = (n - thisDegeneracy == 0) ? new NullSinglet()
                    : new PermutationDegenerateIterator(n - thisDegeneracy);
            subPermutation.setDegeneracy(subDegeneracy);
        }
    }

    public int[] getDegeneracy() {
        return degeneracy.clone();
    }

    /**
     * Readies the iterator for iteration. Must be called after any calls to
     * setDegeneracy.
     */
    public void reset() {
        subPermutation.reset();
        subIndex = subPermutation.next();
        indexIterator.reset();
        if (indexIterator.hasNext()) {
            iLast = indexIterator.next();
            hasNext = true;
        } else {
            hasNext = false;
        }
    }

    public boolean hasNext() {
        return hasNext;
    }

    /**
     * Returns the next iterate, or null if hasNext is false.
     */
    public int[] next() {
        if (!hasNext) {
            return null;
        }
        System.arraycopy(zero, 0, nextIndex, 0, n);
        for (int i = 0; i < iLast.length; i++) {
            nextIndex[iLast[i]] = element;
        }
        int j = 0;
        for (int i = 0; i < n; i++) {
            if (nextIndex[i] == element)
                continue;
            else
                nextIndex[i] = subIndex[j++];
        }
        if (indexIterator.hasNext()) {
            iLast = indexIterator.next();
        } else {
            if (subPermutation.hasNext()) {
                subIndex = subPermutation.next();
                indexIterator.reset();
                iLast = indexIterator.next();
            } else {
                hasNext = false;
            }
        }
        return nextIndex;
    }
    
    public int getD() {
        return n;
    }

    // Iterator that returns one null value (which isn't used) then expires.
    // Used to close recursive implemention of PermutationDegenerateIterator.
    private static class NullSinglet extends PermutationDegenerateIterator {
        private boolean hasNext = false;
        private static final long serialVersionUID = 1L;

        public boolean hasNext() {
            return hasNext;
        }

        public void reset() {
            hasNext = true;
        }

        public int[] next() {
            hasNext = false;
            return null;
        }
    }

    public static void main(String[] args) {
        // int n = Integer.parseInt(args[0]);
        int n = 5;
        PermutationDegenerateIterator p = new PermutationDegenerateIterator(n);
        p.setDegeneracy(new int[] { 3, 1, 1 });
        int count = 0;
        p.reset();
        while (p.hasNext) {
            int[] a = p.next();
            count++;
            System.out.println(count + ". " + Arrays.toString(a));
            // System.out.print("{");
            // for(int i=0; i<a.length-1; i++) System.out.print(a[i]+", ");
            // System.out.println(a[a.length-1]+"}");
        }
        long expectedCount = SpecialFunctions.factorial(n);
        int[] deg = p.getDegeneracy();
        for (int i = 0; i < deg.length; i++) {
            expectedCount /= SpecialFunctions.factorial(deg[i]);
        }
        System.out.println("Expected count: " + expectedCount);
    }

    private boolean hasNext = false;
    private final int n;
    private final int[] nextIndex;
    private final int[] zero;

    // the subPermutation iterator returns permutations that of elements that do
    // not include the one added by this iterator.
    private PermutationDegenerateIterator subPermutation;

    // subIndex is the array returned by the subPermutation iterator. The value
    // "element" is inserted in this
    // array at the positions indicated by iLast to produce the int[] that is
    // returned by this iterator
    private int[] subIndex;

    // a copy of degeneracy is kept only to service the getDegeneracy method
    private int[] degeneracy;

    // this iterator is responsible to placing the value "element" in the
    // returned int[] array
    private int element;

    // indexIterator is used to determine the places where this iterator's
    // element is inserted into the subiterator array
    private IndexIteratorTriangular indexIterator;

    // iLast is generated by indexIterator and points to the places where
    // element should be inserted
    private int[] iLast;

    private static final long serialVersionUID = 1L;

}
