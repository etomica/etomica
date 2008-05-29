package etomica.atom;

import etomica.api.IAtom;
import etomica.api.IAtomSet;

public class AtomArrayListTemp implements IAtomSet {

	private float trimThreshold = 0.8f;
	private IAtom[] atomList;
	private static int DEFAULT_INIT_SIZE = 20;
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
			atomList = null;
			atomList = tempList;
		}
	}

	public void maybeTrimToSize() {
		float usedSpace = (float)itemsInList / (float)atomList.length;
		if(usedSpace < trimThreshold &&
		   itemsInList < atomList.length) {
			IAtom[] tempList = toArray();
			itemsInList = tempList.length;
			atomList = null;
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
		boolean empty = false;
		if(itemsInList == 0) {
			empty = true;
		}
		return empty;
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
		IAtom[] tempList = null;
		if(itemsInList > 0) {
			tempList = new IAtom[itemsInList];
			for(int i = 0; i < itemsInList; i++) {
				tempList[i] = atomList[i];
			}
		}
		return tempList;
	}

	public IAtom set(int index, IAtom element) {
		IAtom oldAtom = null;
		if(index > 0 && index < itemsInList) {
			oldAtom = atomList[index];
			atomList[index] = element;
		}
		else {
			throw new IndexOutOfBoundsException("AtomArrayList.set index out of bounds");
		}
		
		return oldAtom;
	}

	public boolean add(IAtom atom) {
//		true if this collection changed as a result of the call
		if(itemsInList == atomList.length) {
			IAtom[] tempList = new IAtom[atomList.length+1];
			for(int i = 0; i < atomList.length; i++) {
				tempList[i] = atomList[i];
			}
			atomList = null;
			atomList = tempList;
		}
		atomList[itemsInList] = atom; 
		itemsInList++;
		
		return true;
	}

	public void addAll(IAtomSet atoms) {
		if((itemsInList + atoms.getAtomCount()) > atomList.length) {
			IAtom[] tempList = new IAtom[itemsInList + atoms.getAtomCount()];
			for(int i = 0; i < atomList.length; i++) {
				tempList[i] = atomList[i];
			}
			atomList = null;
			atomList = tempList;
		}
		for(int i = 0; i < atoms.getAtomCount(); i++) {
			atomList[itemsInList] = atoms.getAtom(i);
			itemsInList++;
		}

	}

	public IAtom remove(int index) {
		IAtom atom = null;
		
		if(index >= 0 && index < itemsInList) {
			atom = atomList[index];
			for(int i = index; i < itemsInList-1; i++) {
				atomList[i] = atomList[i+1];
			}
			atomList[itemsInList-1] = null;
		    itemsInList--;
		}
		else {
			throw new IndexOutOfBoundsException("AtomArrayList.remove invalid index");
		}
		return atom;
	}

	public IAtom removeAndReplace(int index) {
		IAtom atom = null;
		
		if(index >= 0 && index < itemsInList) {
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
		else {
			throw new IndexOutOfBoundsException("AtomArrayList.remove invalid index");
		}
		return atom;
	}

	public void clear() {
		int size = atomList.length;
		atomList = null;
		atomList = new IAtom[size];
		itemsInList = 0;
	}

	public int sizeOfArray() {
		return atomList.length;
	}

	public int getAtomCount() {
		return itemsInList;
	}
	
	public IAtom getAtom(int index) {
		return atomList[index];
	}

}
