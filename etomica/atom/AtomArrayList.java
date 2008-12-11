package etomica.atom;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;
import etomica.util.Debug;

public class AtomArrayList implements IAtomList {

    protected float trimThreshold = 0.8f;
    protected IAtomLeaf[] atomList;
    protected static int DEFAULT_INIT_SIZE = 20;
    protected static float SIZE_INCREASE_RATIO = 0.3f;
    protected int itemsInList = 0;

    public AtomArrayList() {
        atomList = new IAtomLeaf[DEFAULT_INIT_SIZE];
    }

    public AtomArrayList(int initialSize) {
        atomList = new IAtomLeaf[initialSize];
    }

    public static float getSizeIncreaseRatio() {
        return SIZE_INCREASE_RATIO;
    }

    public void trimToSize() {
        if(itemsInList < atomList.length) {
            IAtomLeaf[] tempList = toArray();
            itemsInList = tempList.length;
            atomList = tempList;
        }
    }

    public void maybeTrimToSize() {
        if(itemsInList < trimThreshold * atomList.length) {
            IAtomLeaf[] tempList = toArray();
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
            IAtomLeaf[] tempList = new IAtomLeaf[minCapacity];
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

    protected IAtomLeaf[] toArray() {
        IAtomLeaf[] tempList = new IAtomLeaf[itemsInList];

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

    public int indexOf(IAtomLeaf elem) {
        for(int i = 0; i < itemsInList; i++) {
            if(elem == atomList[i]) {
                return i;
            }
        }
        return -1;
    }

    public IAtomLeaf[] toAtomLeafArray() {
        IAtomLeaf[] tempList = new IAtomLeaf[itemsInList];

        for(int i = 0; i < itemsInList; i++) {
            tempList[i] = atomList[i];
        }
        return tempList;
    }

    public IAtomLeaf set(int index, IAtomLeaf element) {
        IAtomLeaf oldAtom = null;
        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("AtomLeafArrayList.set index out of bounds");
        }

        oldAtom = atomList[index];
        atomList[index] = element;

        return oldAtom;
    }

    public boolean add(IAtomLeaf atom) {

        if(itemsInList == atomList.length) {
            IAtomLeaf[] tempList = new IAtomLeaf[(int)((float)itemsInList * (1.0f + SIZE_INCREASE_RATIO)+1)];

            for(int i = 0; i < atomList.length; i++) {
                tempList[i] = atomList[i];
            }
            atomList = tempList;
        }
        atomList[itemsInList] = atom; 
        itemsInList++;

        return true;
    }

    public void addAll(IAtomList atoms) {
        if((itemsInList + atoms.getAtomCount()) > atomList.length) {
            IAtomLeaf[] tempList = new IAtomLeaf[(int)((float)itemsInList * (1.0f + SIZE_INCREASE_RATIO)) +
                                         atoms.getAtomCount()];
            for(int i = 0; i < atomList.length; i++) {
                tempList[i] = atomList[i];
            }
            atomList = tempList;
        }
        for(int i = 0; i < atoms.getAtomCount(); i++) {
            atomList[itemsInList] = atoms.getAtom(i);
            itemsInList++;
        }
    }

    public IAtomLeaf remove(int index) {
        IAtomLeaf atom = null;

        if(index < 0 || index >= itemsInList) {
            throw new IndexOutOfBoundsException("AtomLeafArrayList.remove invalid index");
        }

        atom = atomList[index];
        for(int i = index; i < itemsInList-1; i++) {
            atomList[i] = atomList[i+1];
        }
        atomList[itemsInList-1] = null;
        itemsInList--;

        return atom;
    }

    public IAtomLeaf removeAndReplace(int index) {
        IAtomLeaf atom = null;

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
    
    public IAtomLeaf getAtom(int index) {
        if (Debug.ON && (index < 0 || index >= itemsInList)) {
            throw new IndexOutOfBoundsException("AtomLeafArrayList.remove invalid index");
        }
        return atomList[index];
    }
}
