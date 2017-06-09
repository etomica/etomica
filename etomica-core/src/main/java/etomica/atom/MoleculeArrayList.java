/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.util.Debug;

public class MoleculeArrayList implements IMoleculeList {

    protected float trimThreshold = 0.8f;
    protected IMolecule[] molecules;
    protected static int DEFAULT_INIT_SIZE = 20;
    protected static float SIZE_INCREASE_RATIO = 0.3f;
    protected int itemsInList = 0;

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
            IMolecule[] tempList = new IMolecule[minCapacity];
            for(int i = 0; i < itemsInList; i++) {
                tempList[i] = molecules[i];
            }
            molecules = null;
            molecules = tempList;
        }
    }

    public boolean isEmpty() {
        return itemsInList == 0;
    }

    protected IMolecule[] toArray() {
        IMolecule[] tempList = new IMolecule[itemsInList];

        for(int i = 0; i < itemsInList; i++) {
            tempList[i] = molecules[i];
        }
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

        for(int i = 0; i < itemsInList; i++) {
            tempList[i] = molecules[i];
        }
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

    public boolean add(IMolecule atom) {

        if(itemsInList == molecules.length) {
            IMolecule[] tempList = new IMolecule[(int)((float)itemsInList * (1.0f + SIZE_INCREASE_RATIO)+1)];

            for(int i = 0; i < molecules.length; i++) {
                tempList[i] = molecules[i];
            }
            molecules = tempList;
        }
        molecules[itemsInList] = atom; 
        itemsInList++;

        return true;
    }

    public void addAll(IMoleculeList atoms) {
        if((itemsInList + atoms.getMoleculeCount()) > molecules.length) {
            IMolecule[] tempList = new IMolecule[(int)((float)itemsInList * (1.0f + SIZE_INCREASE_RATIO)) +
                                         atoms.getMoleculeCount()];
            for(int i = 0; i < molecules.length; i++) {
                tempList[i] = molecules[i];
            }
            molecules = tempList;
        }
        for(int i = 0; i < atoms.getMoleculeCount(); i++) {
            molecules[itemsInList] = atoms.getMolecule(i);
            itemsInList++;
        }
    }

    public IMolecule remove(int index) {
        IMolecule atom = null;

        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("MoleculeArrayList.remove invalid index");
        }

        atom = molecules[index];
        for(int i = index; i < itemsInList-1; i++) {
            molecules[i] = molecules[i+1];
        }
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

    public int getMoleculeCount() {
        return itemsInList;
    }
    
    public IMolecule getMolecule(int index) {
        if (Debug.ON && (index < 0 || index >= itemsInList)) {
            throw new IndexOutOfBoundsException("MoleculeArrayList.remove invalid index");
        }
        return molecules[index];
    }
}
