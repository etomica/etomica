/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

public class Int1Darray {

    public int offset;
    public final int[] storage;
    public Int1Darray(int[] storage, int offset) {
        this.storage = storage;
        this.offset = offset;
    }

    public void verifyStorage(int[] storage) {
        if (this.storage != storage) throw new RuntimeException("wrong storage");
    }

    public void setOffset(int offset) {
        this.offset = offset;
    }

    public int getOffset() {
        return offset;
    }

    public int get(int i) {
        return storage[offset + i];
    }

    public void set(int i, int x) {
        storage[offset + i] = x;
    }
}
