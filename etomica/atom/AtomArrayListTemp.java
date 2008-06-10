package etomica.atom;

import etomica.api.IAtom;
import etomica.api.IAtomSet;

public class AtomArrayListTemp implements IAtomSet {

	private float trimThreshold = 0.8f;
	private IAtom[] atomList;
	private static int DEFAULT_INIT_SIZE = 20;
	private static float SIZE_INCREASE_RATIO = 0.3f;
	int itemsInList = 0;

	public AtomArrayListTemp() {
		atomList = new IAtom[DEFAULT_INIT_SIZE];
	}

	public AtomArrayListTemp(int initialSize) {
		atomList = new IAtom[initialSize];
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
			IAtom[] tempList = new IAtom[minCapacity];
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

	public int indexOf(IAtom elem) {
		int index = -1;
		for(int i = 0; i < itemsInList; i++) {
			if(elem == atomList[i]) {
				index = i;
				break;
			}
		}
		return index;
	}

	public IAtom[] toArray() {
		IAtom[] tempList = new IAtom[itemsInList];

		for(int i = 0; i < itemsInList; i++) {
			tempList[i] = atomList[i];
		}
		return tempList;
	}

	public IAtom set(int index, IAtom element) {
		IAtom oldAtom = null;
		if(index < 0 || index >= itemsInList) {
			throw new IndexOutOfBoundsException("AtomArrayList.set index out of bounds");
		}
		else {
			oldAtom = atomList[index];
			atomList[index] = element;
		}
		
		return oldAtom;
	}

	public boolean add(IAtom atom) {

		if(itemsInList == atomList.length) {
			IAtom[] tempList = new IAtom[(int)((float)itemsInList * (1.0f + SIZE_INCREASE_RATIO))];

			for(int i = 0; i < atomList.length; i++) {
				tempList[i] = atomList[i];
			}
            atomList = tempList;
		}
		atomList[itemsInList] = atom; 
		itemsInList++;
		
		return true;
	}

	public void addAll(IAtomSet atoms) {
		if((itemsInList + atoms.getAtomCount()) > atomList.length) {
			IAtom[] tempList = new IAtom[(int)((float)itemsInList * (1.0f + SIZE_INCREASE_RATIO)) +
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

	public IAtom remove(int index) {
		IAtom atom = null;
		
		if(index < 0 || index >= itemsInList) {
			throw new IndexOutOfBoundsException("AtomArrayList.remove invalid index");
		}
		else {
			atom = atomList[index];
			for(int i = index; i < itemsInList-1; i++) {
				atomList[i] = atomList[i+1];
			}
			atomList[itemsInList-1] = null;
		    itemsInList--;
		}

		return atom;
	}

	public IAtom removeAndReplace(int index) {
		IAtom atom = null;

		if(index < 0 || index >= itemsInList) {
			throw new IndexOutOfBoundsException("AtomArrayList.remove invalid index");
		}
		else {
			atom = atomList[index];
			if(index == atomList.length-1) {
				atomList[index] = null;
			}
			else {
				atomList[index] = atomList[itemsInList-1];
				atomList[itemsInList-1] = null;
			}
			itemsInList--;
		}

		return atom;
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

	public int getAtomCount() {
		return itemsInList;
	}
	
	public IAtom getAtom(int index) {
		if(index < 0 || index >= itemsInList) {
			throw new IndexOutOfBoundsException("AtomArrayList.remove invalid index");
		}
		return atomList[index];
	}

}
