/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

import java.util.Arrays;

public class Int2Darray {

    public int[] offsets;
    public int[][] storage;
    public Int2Darray() {
        storage = new int[0][];
        offsets = new int[0];
    }

    public void setIndex(int idx, int[] storage, int offset) {
        if (offsets.length <= idx) {
            offsets = Arrays.copyOf(offsets, idx+1);
            this.storage = Arrays.copyOf(this.storage, idx+1);
        }
        this.storage[idx] = storage;
        offsets[idx] = offset;
    }

    public void copyOffset(int idxSrc, int idxDest) {
        storage[idxDest] = storage[idxSrc];
        offsets[idxDest] = offsets[idxSrc];
    }

    public void setOffsets(int[] storage, int[] offsets) {
        this.storage = new int[offsets.length][];
        Arrays.fill(this.storage, storage);
        this.offsets = offsets;
    }

    public int[] getStorage(int idx) {
        return storage[idx];
    }
    public int getOffset(int idx) {
        return offsets[idx];
    }

    public int get(int i, int j) {
        return storage[i][offsets[i] + j];
    }

    public void set(int i, int j, int x) {
        storage[i][offsets[i] + j] = x;
    }
}
