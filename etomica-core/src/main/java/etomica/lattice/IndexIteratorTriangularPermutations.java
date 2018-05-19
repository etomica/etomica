/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

import etomica.math.discrete.PermutationDegenerateIterator;
import java.util.Arrays;

/**
 * Generates all unique iterates formed as permutations of iterates given by IndexIteratorTriangular with
 * allowEqualElements set to true.   This aims to generate iterates over a pyrimidal region, such that each subsequent
 * iterate has a sum-squared magnitude that is greater than that of previous iterates (this is accomplished
 * only approximately; some iterates are "closer" to origin than some others that precede it).  This ordering
 * contrasts with IndexIteratorRectangular which nests iteration over each coordinate direction.
 * 
 * For example, with D set to 3 (at construction), and maxElement = 2, the following iterates are generated:
 * {0,0,0}<br>
 * {0,0,1}<br>
 * {0,1,0}<br>
 * {1,0,0}<br>
 * {0,1,1}<br>
 * {1,0,1}<br>
 * {1,1,0}<br>
 * {1,1,1}<br>
 * {0,0,2}<br>
 * {0,2,0}<br>
 * {2,0,0}<br>
 * {0,1,2}<br>
 * {1,0,2}<br>
 * {1,2,0}<br>
 * {0,2,1}<br>
 * {2,0,1}<br>
 * {2,1,0}<br>
 * {1,1,2}<br>
 * {1,2,1}<br>
 * {2,1,1}<br>
 * {0,2,2}<br>
 * {2,0,2}<br>
 * {2,2,0}<br>
 * {1,2,2}<br>
 * {2,1,2}<br>
 * {2,2,1}<br>
 * {2,2,2}<br>
 * 
 * @author David Kofke
 *
 */
public class IndexIteratorTriangularPermutations implements IndexIterator, java.io.Serializable {

    public IndexIteratorTriangularPermutations(int D) {
        if(D < 0) throw new IllegalArgumentException("Iterator must have non-negative dimension. Given value is "+D);
        this.D = D;
        coreIterator = new IndexIteratorTriangular(D);
        coreIterator.setAllowEqualElements(true);
        permutationIterator = new PermutationDegenerateIterator(D);
        nextIterate = new int[D];
        vals = new int[D];
        counts = new int[D];
        
    }
    public int getD() {
        return coreIterator.getD();
    }

    public boolean hasNext() {
        return hasNext;
    }

    public int[] next() {
        for(int i=0; i<nextIterate.length; i++) {
            nextIterate[i] = vals[permute[i]];
        }
        if(permutationIterator.hasNext()) {
            permute = permutationIterator.next();
        } else {
            if(coreIterator.hasNext()) {
                resetPermutationIterator(coreIterator.next());
                permute = permutationIterator.next();
            } else {
                hasNext = false;
            }
        }
        return nextIterate;
    }

    public void reset() {
        hasNext = false;
        coreIterator.reset();
        if(coreIterator.hasNext()) {
            resetPermutationIterator(coreIterator.next());
            if(permutationIterator.hasNext()) {
                permute = permutationIterator.next();
                hasNext = true;
            }
        }
    }
    
    public IndexIteratorTriangular getCoreIterator() {
        return coreIterator;
    }
    
    private void resetPermutationIterator(int[] baseIterate) {
        vals[0] = baseIterate[0];
        counts[0] = 1;
        int j = 0;
        for(int i=1; i<D; i++) {
            if(baseIterate[i] == vals[j]) {
                counts[j]++;
            } else {
                j++;
                vals[j] = baseIterate[i];
                counts[j] = 1;
            }
        }
        permutationIterator.setDegeneracy(counts, j+1);
        permutationIterator.reset();
    }
    
    /**
     * Method to test and demonstrate class.
     */
    public static void main(String[] args) {
        IndexIteratorTriangularPermutations iterator = new IndexIteratorTriangularPermutations(1);
        iterator.getCoreIterator().setMaxElement(2);
        iterator.reset();
        System.out.println("Start");
        int count = 0;
        while(iterator.hasNext()) {
            int[] a = iterator.next();
            int sum = 0;
            for(int i=0; i<a.length; i++) {
                sum += a[i]*a[i];
            }
            System.out.println(++count+". "+Arrays.toString(a) +" "/*"<br>");/*/+Math.sqrt(sum));
        }
    }

    private final IndexIteratorTriangular coreIterator;
    private final PermutationDegenerateIterator permutationIterator;
    private int[] permute;
    private int[] nextIterate;
    private final int[] vals;
    private final int[] counts;
    private final int D;
    private boolean hasNext;
    
    private static final long serialVersionUID = 1L;
}
