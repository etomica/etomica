/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 *
 * Copyright (C) 2005  Miguel, Jmol Development, www.jmol.org
 *
 * Contact: jmol-developers@lists.sf.net
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

package etomica.util.collections;

import java.util.Arrays;

/**
 * This class adds clear and remove methods to jmol's implementation.
 * <p>
 * The put method attempts grows the entries array if the average number of
 * entries per hash exceeds 1.
 */
public class Int2IntHash {
    protected int entryCount;
    protected Entry[] entries;

    public Int2IntHash(int initialCapacity) {
        entries = new Entry[initialCapacity];
    }

    public void clear() {
        Arrays.fill(entries, null);
    }

    public void remove(int key) {
        Entry[] entries = this.entries;
        int hash = (key & 0x7FFFFFFF) % entries.length;
        Entry ePrev = null;
        for (Entry e = entries[hash]; e != null; e = e.next) {
            if (e.key == key) {
                if (ePrev == null) {
                    entries[hash] = e.next;
                } else {
                    ePrev.next = e.next;
                }
                return;
            }
            ePrev = e;
        }
    }

    public int get(int key) {
        Entry[] entries = this.entries;
        int hash = (key & 0x7FFFFFFF) % entries.length;
        for (Entry e = entries[hash]; e != null; e = e.next)
            if (e.key == key)
                return e.value;
        return Integer.MIN_VALUE;
    }

    public void put(int key, int value) {
        Entry[] entries = this.entries;
        int hash = (key & 0x7FFFFFFF) % entries.length;
        for (Entry e = entries[hash]; e != null; e = e.next)
            if (e.key == key) {
                e.value = value;
                return;
            }
        if (entryCount > entries.length)
            rehash();
        entries = this.entries;
        hash = (key & 0x7FFFFFFF) % entries.length;
        entries[hash] = new Entry(key, value, entries[hash]);
        ++entryCount;
    }

    private void rehash() {
        Entry[] oldEntries = entries;
        int oldSize = oldEntries.length;
        int newSize = oldSize * 2 + 1;
        Entry[] newEntries = new Entry[newSize];

        for (int i = oldSize; --i >= 0; ) {
            for (Entry e = oldEntries[i]; e != null; ) {
                Entry t = e;
                e = e.next;

                int hash = (t.key & 0x7FFFFFFF) % newSize;
                t.next = newEntries[hash];
                newEntries[hash] = t;
            }
        }
        entries = newEntries;
    }

    static class Entry {
        int key;
        int value;
        Entry next;

        Entry(int key, int value, Entry next) {
            this.key = key;
            this.value = value;
            this.next = next;
        }
    }
}