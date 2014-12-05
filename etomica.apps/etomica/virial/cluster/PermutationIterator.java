/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import java.util.Arrays;

public class PermutationIterator {

    public PermutationIterator(int[] indices) {
        this.indices = indices;
        Arrays.sort(this.indices);
        nIndices = indices.length;
        rIndices = new int[nIndices];
        int nUniqueIndices = 1;
        for (int i=1; i<indices.length; i++) {
            if (indices[i] != indices[i-1]) {
                nUniqueIndices++;
            }
        }
        maxIndex = indices[indices.length-1];
        nEachIndices = new int[maxIndex+1];
        for (int i=0; i<nIndices; i++) {
            nEachIndices[indices[i]]++;
        }
        nEachIndicesCheck = new int[maxIndex+1];
    }
    
    public void reset() {
        for (int i=0; i<nIndices; i++) {
            rIndices[i] = 0;
        }
        rIndices[nIndices-1] = -1;
    }
    
    public int[] next() {
        int i = nIndices-1;
        while (true) {
            rIndices[i]++;
            boolean flag = false;
            while (rIndices[i] > maxIndex) {
                flag = true;
                if (i == 0) {
                    return null;
                }
                rIndices[i] = 0;
                i--;
                rIndices[i]++;
            }
            if (flag) {
                i = nIndices-1;
            }
            // now check # of each index
            for (int j=0; j<nEachIndicesCheck.length; j++) {
                nEachIndicesCheck[j] = 0;
            }
            for (int j=0; j<nIndices; j++) {
                nEachIndicesCheck[rIndices[j]]++;
            }
            boolean success = true;
            for (int j=0; j<nEachIndicesCheck.length; j++) {
                if (nEachIndicesCheck[j] != nEachIndices[j]) {
                    success = false;
                    break;
                }
            }
            if (success) break;
        }
        return rIndices;
    }
    
    protected final int[] indices;
    protected final int[] rIndices;
    protected final int maxIndex;
    protected final int nIndices;
    protected final int[] nEachIndices, nEachIndicesCheck;
    
    public static void main(String[] args) {
        PermutationIterator it = new PermutationIterator(new int[]{0,0,1,1});
        it.reset();
        while (true) {
            int[] i = it.next();
            if (i == null) break;
            System.out.println(Arrays.toString(i));
        }
    }
}
