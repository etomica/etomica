/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.util.Debug;

import java.util.Arrays;

public final class AtomArrayList implements IAtomList {

    private float trimThreshold = 0.8f;
    private IAtom[] atomList;
    private static final int DEFAULT_INIT_SIZE = 20;
    private static final float SIZE_INCREASE_RATIO = 0.3f;
    private int itemsInList = 0;

    public AtomArrayList() {
        atomList = new IAtom[DEFAULT_INIT_SIZE];
    }

    public AtomArrayList(int initialSize) {
        atomList = new IAtom[initialSize];
    }

    public static float getSizeIncreaseRatio() {
        return SIZE_INCREASE_RATIO;
    }

    public void trimToSize() {
        if(itemsInList < atomList.length) {
            IAtom[] tempList = toArray();
            itemsInList = tempList.length;
            atomList = tempList;
        }
    }

    public void maybeTrimToSize() {
        if(itemsInList < trimThreshold * atomList.length) {
            IAtom[] tempList = toArray();
            itemsInList = tempList.length;
            atomList = tempList;
        }
    }

    public float getTrimThreshold() {
        return trimThreshold;
    }

    public void setTrimThreshold(float newTrimThreshold) {
        trimThreshold = newTrimThreshold;
    }

    public void ensureCapacity(int minCapacity) {
        if(minCapacity > atomList.length) {
            atomList = Arrays.copyOf(atomList, minCapacity);
        }
    }

    public boolean isEmpty() {
        return itemsInList == 0;
    }

    protected IAtom[] toArray() {
        IAtom[] tempList = new IAtom[itemsInList];

        System.arraycopy(atomList, 0, tempList, 0, itemsInList);
        return tempList;
    }

    public void clear() {
        for(int i = 0; i < itemsInList; i++) {
            atomList[i] = null;
        }
        itemsInList = 0;
    }

    public int sizeOfArray() {
        return atomList.length;
    }

    public int indexOf(IAtom elem) {
        for(int i = 0; i < itemsInList; i++) {
            if(elem == atomList[i]) {
                return i;
            }
        }
        return -1;
    }
    
    /**
     * Returns true if the given atom is in the list, false otherwise
     */
    public boolean contains(IAtom elem) {
        return (indexOf(elem) != -1);
    }

    public IAtom[] toAtomLeafArray() {
        IAtom[] tempList = new IAtom[itemsInList];

        System.arraycopy(atomList, 0, tempList, 0, itemsInList);
        return tempList;
    }

    public IAtom set(int index, IAtom element) {
        IAtom oldAtom = null;
        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("AtomLeafArrayList.set index out of bounds");
        }

        oldAtom = atomList[index];
        atomList[index] = element;

        return oldAtom;
    }

    public boolean add(IAtom atom) {
        if(itemsInList == atomList.length) {
            atomList = Arrays.copyOf(atomList, (int)((float)itemsInList * (1.0f + SIZE_INCREASE_RATIO)+1));
        }
        atomList[itemsInList] = atom; 
        itemsInList++;

        return true;
    }

    public void addAll(IAtomList atoms) {
        if((itemsInList + atoms.getAtomCount()) > atomList.length) {
            atomList = Arrays.copyOf(
                    atomList,
                    (int)((float)itemsInList * (1.0f + SIZE_INCREASE_RATIO)) + atoms.getAtomCount()
            );
        }
        if (atoms instanceof AtomArrayList) {
            System.arraycopy(((AtomArrayList) atoms).atomList, 0, this.atomList, itemsInList, atoms.getAtomCount());
            itemsInList += atoms.getAtomCount();
        } else {
            for(int i = 0; i < atoms.getAtomCount(); i++) {
                atomList[itemsInList] = atoms.getAtom(i);
                itemsInList++;
            }
        }
    }

    /**
     * This removes the atom at the given index and shifts the remaining atoms
     * back (maintaining the order).
     */
    public IAtom remove(int index) {
        IAtom atom = null;

        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("AtomLeafArrayList.remove invalid index");
        }

        atom = atomList[index];
        System.arraycopy(atomList, index + 1, atomList, index, itemsInList - 1 - index);
        atomList[itemsInList-1] = null;
        itemsInList--;

        return atom;
    }

    /**
     * This removes the atom at the given index and replaces the hole with the
     * last atom.  This is faster than the remove() method.
     */
    public IAtom removeAndReplace(int index) {
        IAtom atom = null;

        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("AtomLeafArrayList.remove invalid index");
        }

        atom = atomList[index];
        if(index == atomList.length-1) {
            atomList[index] = null;
        }
        else {
            atomList[index] = atomList[itemsInList-1];
            atomList[itemsInList-1] = null;
        }
        itemsInList--;

        return atom;
    }

    public int getAtomCount() {
        return itemsInList;
    }
    
    public IAtom getAtom(int index) {
        if (Debug.ON && (index < 0 || index >= itemsInList)) {
            throw new IndexOutOfBoundsException("AtomLeafArrayList.remove invalid index");
        }
        return atomList[index];
    }

    public String toString() {
        if (itemsInList == 0) return "";
        String str = ""+atomList[0];
        for (int i=1; i<itemsInList; i++) {
            str += " "+atomList[i];
        }
        return str;
    }
}
