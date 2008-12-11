package etomica.atom;

import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.util.Debug;

public class MoleculeArrayList implements IMoleculeList {

    protected float trimThreshold = 0.8f;
    protected IMolecule[] atomList;
    protected static int DEFAULT_INIT_SIZE = 20;
    protected static float SIZE_INCREASE_RATIO = 0.3f;
    protected int itemsInList = 0;

    public MoleculeArrayList() {
        atomList = new IMolecule[DEFAULT_INIT_SIZE];
    }

    public MoleculeArrayList(int initialSize) {
        atomList = new IMolecule[initialSize];
    }

    public static float getSizeIncreaseRatio() {
        return SIZE_INCREASE_RATIO;
    }

    public void trimToSize() {
        if(itemsInList < atomList.length) {
            IMolecule[] tempList = toArray();
            itemsInList = tempList.length;
            atomList = tempList;
        }
    }

    public void maybeTrimToSize() {
        if(itemsInList < trimThreshold * atomList.length) {
            IMolecule[] tempList = toArray();
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
            IMolecule[] tempList = new IMolecule[minCapacity];
            for(int i = 0; i < itemsInList; i++) {
                tempList[i] = atomList[i];
            }
            atomList = null;
            atomList = tempList;
        }
    }

    public boolean isEmpty() {
        return itemsInList == 0;
    }

    protected IMolecule[] toArray() {
        IMolecule[] tempList = new IMolecule[itemsInList];

        for(int i = 0; i < itemsInList; i++) {
            tempList[i] = atomList[i];
        }
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

    public int indexOf(IMolecule elem) {
        for(int i = 0; i < itemsInList; i++) {
            if(elem == atomList[i]) {
                return i;
            }
        }
        return -1;
    }

    public IMolecule[] toMoleculeArray() {
        IMolecule[] tempList = new IMolecule[itemsInList];

        for(int i = 0; i < itemsInList; i++) {
            tempList[i] = atomList[i];
        }
        return tempList;
    }

    public IMolecule set(int index, IMolecule element) {
        IMolecule oldAtom = null;
        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("MoleculeArrayList.set index out of bounds");
        }

        oldAtom = atomList[index];
        atomList[index] = element;

        return oldAtom;
    }

    public boolean add(IMolecule atom) {

        if(itemsInList == atomList.length) {
            IMolecule[] tempList = new IMolecule[(int)((float)itemsInList * (1.0f + SIZE_INCREASE_RATIO)+1)];

            for(int i = 0; i < atomList.length; i++) {
                tempList[i] = atomList[i];
            }
            atomList = tempList;
        }
        atomList[itemsInList] = atom; 
        itemsInList++;

        return true;
    }

    public void addAll(IMoleculeList atoms) {
        if((itemsInList + atoms.getMoleculeCount()) > atomList.length) {
            IMolecule[] tempList = new IMolecule[(int)((float)itemsInList * (1.0f + SIZE_INCREASE_RATIO)) +
                                         atoms.getMoleculeCount()];
            for(int i = 0; i < atomList.length; i++) {
                tempList[i] = atomList[i];
            }
            atomList = tempList;
        }
        for(int i = 0; i < atoms.getMoleculeCount(); i++) {
            atomList[itemsInList] = atoms.getMolecule(i);
            itemsInList++;
        }
    }

    public IMolecule remove(int index) {
        IMolecule atom = null;

        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("MoleculeArrayList.remove invalid index");
        }

        atom = atomList[index];
        for(int i = index; i < itemsInList-1; i++) {
            atomList[i] = atomList[i+1];
        }
        atomList[itemsInList-1] = null;
        itemsInList--;

        return atom;
    }

    public IMolecule removeAndReplace(int index) {
        IMolecule atom = null;

        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("MoleculeArrayList.remove invalid index");
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

    public int getMoleculeCount() {
        return itemsInList;
    }
    
    public IMolecule getMolecule(int index) {
        if (Debug.ON && (index < 0 || index >= itemsInList)) {
            throw new IndexOutOfBoundsException("MoleculeArrayList.remove invalid index");
        }
        return atomList[index];
    }
}
