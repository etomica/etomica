/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.util.Debug;

import java.util.Arrays;

public final class MoleculeArrayList implements IMoleculeList {

    private float trimThreshold = 0.8f;
    private IMolecule[] molecules;
    private static final int DEFAULT_INIT_SIZE = 20;
    private static final float SIZE_INCREASE_RATIO = 0.3f;
    private int itemsInList = 0;

    public MoleculeArrayList() {
        molecules = new IMolecule[DEFAULT_INIT_SIZE];
    }

    public MoleculeArrayList(int initialSize) {
        molecules = new IMolecule[initialSize];
    }

    public static float getSizeIncreaseRatio() {
        return SIZE_INCREASE_RATIO;
    }

    public void trimToSize() {
        if(itemsInList < molecules.length) {
            IMolecule[] tempList = toArray();
            itemsInList = tempList.length;
            molecules = tempList;
        }
    }

    public void maybeTrimToSize() {
        if(itemsInList < trimThreshold * molecules.length) {
            IMolecule[] tempList = toArray();
            itemsInList = tempList.length;
            molecules = tempList;
        }
    }

    public float getTrimThreshold() {
        return trimThreshold;
    }

    public void setTrimThreshold(float newTrimThreshold) {
        trimThreshold = newTrimThreshold;
    }

    public void ensureCapacity(int minCapacity) {
        if(minCapacity > molecules.length) {
            molecules = Arrays.copyOf(molecules, minCapacity);
        }
    }

    public boolean isEmpty() {
        return itemsInList == 0;
    }

    protected IMolecule[] toArray() {
        IMolecule[] tempList = new IMolecule[itemsInList];

        System.arraycopy(molecules, 0, tempList, 0, itemsInList);
        return tempList;
    }

    public void clear() {
        for(int i = 0; i < itemsInList; i++) {
            molecules[i] = null;
        }
        itemsInList = 0;
    }

    public int sizeOfArray() {
        return molecules.length;
    }

    public int indexOf(IMolecule elem) {
        for(int i = 0; i < itemsInList; i++) {
            if(elem == molecules[i]) {
                return i;
            }
        }
        return -1;
    }

    public IMolecule[] toMoleculeArray() {
        IMolecule[] tempList = new IMolecule[itemsInList];

        System.arraycopy(molecules, 0, tempList, 0, itemsInList);
        return tempList;
    }

    public IMolecule set(int index, IMolecule element) {
        IMolecule oldAtom = null;
        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("MoleculeArrayList.set index out of bounds");
        }

        oldAtom = molecules[index];
        molecules[index] = element;

        return oldAtom;
    }

    public boolean add(IMolecule molecule) {

        if(itemsInList == molecules.length) {
            molecules = Arrays.copyOf(molecules, (int)((float)itemsInList * (1.0f + SIZE_INCREASE_RATIO)+1));
        }
        molecules[itemsInList] = molecule;
        itemsInList++;

        return true;
    }

    public void addAll(IMoleculeList moleculeList) {
        if ((itemsInList + moleculeList.size()) > molecules.length) {
            molecules = Arrays.copyOf(
                    molecules,
                    (int) ((float) itemsInList * (1.0f + SIZE_INCREASE_RATIO)) + moleculeList.size()
            );
        }
        if (moleculeList instanceof MoleculeArrayList) {
            System.arraycopy(((MoleculeArrayList) moleculeList).molecules, 0, this.molecules, itemsInList, moleculeList.size());
            itemsInList += moleculeList.size();
        } else {
            for (int i = 0; i < moleculeList.size(); i++) {
                molecules[itemsInList] = moleculeList.get(i);
                itemsInList++;
            }
        }
    }

    public IMolecule remove(int index) {
        IMolecule atom = null;

        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("MoleculeArrayList.remove invalid index");
        }

        atom = molecules[index];
        System.arraycopy(molecules, index + 1, molecules, index, itemsInList - 1 - index);
        molecules[itemsInList-1] = null;
        itemsInList--;

        return atom;
    }

    public IMolecule removeAndReplace(int index) {
        IMolecule atom = null;

        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("MoleculeArrayList.remove invalid index");
        }

        atom = molecules[index];
        if(index == molecules.length-1) {
            molecules[index] = null;
        }
        else {
            molecules[index] = molecules[itemsInList-1];
            molecules[itemsInList-1] = null;
        }
        itemsInList--;

        return atom;
    }

    public int size() {
        return itemsInList;
    }
    
    public IMolecule get(int index) {
        if (Debug.ON && (index < 0 || index >= itemsInList)) {
            throw new IndexOutOfBoundsException("MoleculeArrayList.remove invalid index");
        }
        return molecules[index];
    }
}
