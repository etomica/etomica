/*
 * @(#)ArrayList.java	1.17 98/09/30
 *
 * Copyright 1997, 1998 by Sun Microsystems, Inc.,
 * 901 San Antonio Road, Palo Alto, California, 94303, U.S.A.
 * All rights reserved.
 *
 * This software is the confidential and proprietary information
 * of Sun Microsystems, Inc. ("Confidential Information").  You
 * shall not disclose such Confidential Information and shall use
 * it only in accordance with the terms of the license agreement
 * you entered into with Sun.
 */

package etomica.atom;

import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.util.Debug;

/**
 * Resizable-array implementation of the List interface.  Implements
 * all optional list operations, and permits all elements, including
 * null.  In addition to implementing the List interface,
 * this class provides methods to manipulate the size of the array that is
 * used internally to store the list.  (This class is roughly equivalent to
 * Vector, except that it is unsynchronized.)
 *
 * The size, isEmpty, get, set,
 * iterator, and listIterator operations run in constant
 * time.  The add operation runs in amortized constant time,
 * that is, adding n elements requires O(n) time.  All of the other operations
 * run in linear time (roughly speaking).  The constant factor is low compared
 * to that for the LinkedList implementation.
 *
 * Each ArrayList instance has a capacity.  The capacity is
 * the size of the array used to store the elements in the list.  It is always
 * at least as large as the list size.  As elements are added an ArrayList,
 * its capacity grows automatically.  The details of the growth policy are not
 * specified beyond the fact that adding an element has constant amortized
 * time cost. 
 *
 * An application can increase the capacity of an ArrayList instance
 * before adding a large number of elements using the ensureCapacity
 * operation.  This may reduce the amount of incremental reallocation.
 *
 * Note that this implementation is not synchronized. If
 * multiple threads access an ArrayList instance concurrently, and at
 * least one of the threads modifies the list structurally, it must be
 * synchronized externally.  (A structural modification is any operation that
 * adds or deletes one or more elements, or explicitly resizes the backing
 * array; merely setting the value of an element is not a structural
 * modification.)  This is typically accomplished by synchronizing on some
 * object that naturally encapsulates the list.  If no such object exists, the
 * list should be "wrapped" using the Collections.synchronizedList
 * method.  This is best done at creation time, to prevent accidental
 * unsynchronized access to the list:
 * 
 *	List list = Collections.synchronizedList(new ArrayList(...));
 * 
 *
 * The iterators returned by this class's iterator and
 * listIterator methods are fail-fast: if list is structurally
 * modified at any time after the iterator is created, in any way except
 * through the iterator's own remove or add methods, the iterator will throw a
 * ConcurrentModificationException.  Thus, in the face of concurrent
 * modification, the iterator fails quickly and cleanly, rather than risking
 * arbitrary, non-deterministic behavior at an undetermined time in the
 * future.
 *
 * @author  Josh Bloch
 * @version 1.17 09/30/98
 * @see	    Collection
 * @see	    List
 * @see	    LinkedList
 * @see	    Vector
 * @see	    Collections#synchronizedList(List)
 * @since JDK1.2
 */

public class AtomArrayList implements Cloneable,
					            java.io.Serializable {
    /**
     * The array buffer into which the elements of the ArrayList are stored.
     * The capacity of the ArrayList is the length of this array buffer.
     */
    private transient Atom elementData[];

    /**
     * The size of the ArrayList (the number of elements it contains).
     *
     * @serial
     */
    private int size;
    
    private float trimThreshold;

    /**
     * Constructs an empty list with the specified initial capacity.
     *
     * @param   initialCapacity   the initial capacity of the list.
     */
    public AtomArrayList(int initialCapacity) {
    	super();
		this.elementData = new Atom[initialCapacity];
        trimThreshold = 0.8F;
    }

    /**
     * Constructs an empty list.
     */
    public AtomArrayList() {
    	this(10);
    }

    /**
     * Trims the capacity of this ArrayList instance to be the
     * list's current size.  An application can use this operation to minimize
     * the storage of an ArrayList instance.
     */
    public void trimToSize() {
		if (size < elementData.length) {
	    	Atom oldData[] = elementData;
	    	elementData = new Atom[size];
	    	System.arraycopy(oldData, 0, elementData, 0, size);
		}
    }
    
    /**
     * Trims the capacity of this AtomArrayList instance (calling 
     * trimToSize) if the fraction of the actual usage is less than 
     * trimThreshold.
     */
    public void maybeTrimToSize() {
        if (size < elementData.length * trimThreshold) {
            trimToSize();
        }
    }

    /**
     * Returns the trim threshhold, the minimum fraction of the array 
     * element usage, below which the array is reallocated in 
     * maybeTrimToSize.  trim threshold defaults to 0.8. 
     */
    public float getTrimThreshold() {
        return trimThreshold;
    }

    /**
     * Sets the trim threshhold, the minimum fraction of the array element 
     * usage, below which the array is reallocated in maybeTrimToSize.
     * trim threshold defaults to 0.8. 
     */
    public void setTrimThreshold(float newTrimThreshold) {
        trimThreshold = newTrimThreshold;
    }

    /**
     * Increases the capacity of this ArrayList instance, if
     * necessary, to ensure  that it can hold at least the number of elements
     * specified by the minimum capacity argument. 
     *
     * @param   minCapacity   the desired minimum capacity.
     */
    public void ensureCapacity(int minCapacity) {
		int oldCapacity = elementData.length;
		if (minCapacity > oldCapacity) {
	    	Atom oldData[] = elementData;
	    	int newCapacity = (oldCapacity * 3)/2 + 1;
    	    if (newCapacity < minCapacity)
    	    	newCapacity = minCapacity;
    	    elementData = new Atom[newCapacity];
	    	System.arraycopy(oldData, 0, elementData, 0, size);
		}
    }

    /**
     * Returns the number of elements in this list.
     *
     * @return  the number of elements in this list.
     */
    public int size() {
    	return size;
    }

    /**
     * Tests if this list has no elements.
     *
     * @return  true if this list has no elements;
     *          false otherwise.
     */
    public boolean isEmpty() {
    	return size == 0;
    }

    /**
     * Returns true if this list contains the specified element.
     *
     * @param o element whose presence in this List is to be tested.
     */
    public boolean contains(Atom elem) {
    	return indexOf(elem) >= 0;
    }

    /**
     * Searches for the first occurence of the given argument, testing 
     * for equality using the equals method. 
     *
     * @param   elem   an atom.
     * @return  the index of the first occurrence of the argument in this
     *          list; returns -1 if the atom is not found.
     * @see     Atom#equals(Atom)
     */
    public int indexOf(Atom elem) {
    	if (elem == null) {
    		for (int i = 0; i < size; i++)
    			if (elementData[i]==null)
    				return i;
    	} else {
    		for (int i = 0; i < size; i++)
    			if (elem.equals(elementData[i]))
    				return i;
    	}
    	return -1;
    }

    /**
     * Returns the index of the last occurrence of the specified atom in
     * this list.
     *
     * @param   elem   the desired element.
     * @return  the index of the last occurrence of the specified atom in
     *          this list; returns -1 if the atom is not found.
     */
    public int lastIndexOf(Atom elem) {
    	if (elem == null) {
    		for (int i = size-1; i >= 0; i--)
    			if (elementData[i]==null)
    				return i;
    	} else {
    		for (int i = size-1; i >= 0; i--)
    			if (elem.equals(elementData[i]))
    				return i;
    	}
    	return -1;
    }

    /**
     * Returns a shallow copy of this ArrayList instance.  (The
     * elements themselves are not copied.)
     *
     * @return  a clone of this ArrayList instance.
     */
    public Object clone() {
    	try { 
    		AtomArrayList v = (AtomArrayList)super.clone();
    		v.elementData = new Atom[size];
    		System.arraycopy(elementData, 0, v.elementData, 0, size);
//	    		v.modCount = 0;
    		return v;
    	} catch (CloneNotSupportedException e) { 
    		// this shouldn't happen, since we are Cloneable
    		throw new InternalError();
    	}
    }

    /**
     * Returns an array containing all of the elements in this list
     * in the correct order.
     *
     * @return an array containing all of the elements in this list
     * 	       in the correct order.
     */
    public Atom[] toArray() {
    	Atom[] result = new Atom[size];
    	System.arraycopy(elementData, 0, result, 0, size);
    	return result;
    }

    /**
     * Returns an array containing all of the elements in this list in the
     * correct order.  The runtime type of the returned array is that of the
     * specified array.  If the list fits in the specified array, it is
     * returned therein.  Otherwise, a new array is allocated with the runtime
     * type of the specified array and the size of this list.
     *
     * If the list fits in the specified array with room to spare (i.e., the
     * array has more elements than the list), the element in the array
     * immediately following the end of the collection is set to
     * null.  This is useful in determining the length of the list
     * only if the caller knows that the list does not contain any
     * null elements.
     *
     * @param a the array into which the elements of the list are to
     *		be stored, if it is big enough; otherwise, a new array of the
     * 		same runtime type is allocated for this purpose.
     * @return an array containing the elements of the list.
     * @throws ArrayStoreException if the runtime type of a is not a supertype
     *         of the runtime type of every element in this list.
     */
    public Atom[] toArray(Atom a[]) {
        if (a.length < size)
        	a = new Atom[size];

        System.arraycopy(elementData, 0, a, 0, size);

        if (a.length > size)
            a[size] = null;

        return a;
    }

    // Positional Access Operations

    /**
     * Returns the element at the specified position in this list.
     * Behavior is undefined if the given index is negative or is larger than the 
     * number of elements in the list.  If Debug.ON is true, an exception is thrown 
     * if the index is out of range.  First index is 0.
     *
     * @param  index index of element to return.
     * @return the element at the specified position in this list.
     */
    public Atom get(int index) {
        RangeCheck(index);
        return elementData[index];
    }

    /**
     * Replaces the element at the specified position in this list with
     * the specified element.
     *
     * @param index index of element to replace.
     * @param element element to be stored at the specified position.
     * @return the element previously at the specified position.
     * @throws    IndexOutOfBoundsException if index out of range
     *		  (index < 0 || index >= size()).
     */
    public Atom set(int index, Atom element) {
    	RangeCheck(index);

    	Atom oldValue = elementData[index];
    	elementData[index] = element;
    	return oldValue;
    }

    /**
     * Appends the specified element to the end of this list.
     *
     * @param o element to be appended to this list.
     * @return true (as per the general contract of Collection.add).
     */
    public boolean add(Atom atom) {
    	ensureCapacity(size + 1);
    	elementData[size++] = atom;
    	return true;
    }

    /**
     * Inserts the specified element at the specified position in this
     * list. Shifts the element currently at that position (if any) and
     * any subsequent elements to the right (adds one to their indices).
     *
     * @param index index at which the specified element is to be inserted.
     * @param element element to be inserted.
     * @throws    IndexOutOfBoundsException if index is out of range
     *		  (index < 0 || index > size()).
     */
    public void add(int index, Atom element) {
    	if (index > size || index < 0)
    		throw new IndexOutOfBoundsException(
    				"Index: "+index+", Size: "+size);

    	ensureCapacity(size+1);
    	System.arraycopy(elementData, index, elementData, index + 1,
    			size - index);
    	elementData[index] = element;
    	size++;
    }
    
    public void addAll(AtomArrayList atoms) {
        ensureCapacity(size+atoms.size());
        int newSize = size + atoms.size();
        for (int i=size; i<newSize; i++) {
            elementData[i] = atoms.get(i-size);
        }
        size = newSize;
    }

    /**
     * Removes the element at the specified position in this list.
     * Shifts any subsequent elements to the left (subtracts one from their
     * indices).
     *
     * @param index the index of the element to be removed.
     * @return the element that was removed from the list.
     * @throws    IndexOutOfBoundsException if index out of range (index
     * 		  < 0 || index >= size()).
     */
    public Atom remove(int index) {
    	RangeCheck(index);

    	Atom oldValue = elementData[index];

    	int numMoved = size - index - 1;
    	if (numMoved > 0)
    		System.arraycopy(elementData, index+1, elementData, index,
    		        numMoved);
    	elementData[--size] = null; // Let gc do its work

    	return oldValue;
    }
    
    /**
     * Removes the element at the specified position in the list.
     * If the element is not the last item in the list, the element is
     * replaced with the last element. 
     * @param index the index of the element to be removed.
     * @return the element that was removed from the list.
     * @throws IndexOutOfBoundsException if the index is out of range
     *         (index < 0 || index >= size()).
     */
    public Atom removeAndReplace(int index) {
        RangeCheck(index);
        
        Atom oldAtom = elementData[index];
        
        size--;
        if (index < size) {
            elementData[index] = elementData[size];
            elementData[size] = null;
        }
        
        return oldAtom;
    }

    /**
     * Removes all of the elements from this list.  The list will
     * be empty after this call returns.
     */
    public void clear() {

    	// Let gc do its work
    	//XXX this is extra work unless the atom is eventually deleted from the system
    	for (int i = 0; i < size; i++)
    		elementData[i] = null;

    	size = 0;
    }

    /**
     * Removes from this List all of the elements whose index is between
     * fromIndex, inclusive and toIndex, exclusive.  Shifts any succeeding
     * elements to the left (reduces their index).
     * This call shortens the list by (toIndex - fromIndex) elements.
     * (If toIndex==fromIndex, this operation has no effect.)
     *
     * @param fromIndex index of first element to be removed.
     * @param fromIndex index after last element to be removed.
     */
    protected void removeRange(int fromIndex, int toIndex) {
    	int numMoved = size - toIndex;
        System.arraycopy(elementData, toIndex, elementData, fromIndex,
        		numMoved);

        // Let gc do its work
        int newSize = size - (toIndex-fromIndex);
        while (size != newSize)
        	elementData[--size] = null;
    }

    /**
     * Check if the given index is in range.  If not, throw an appropriate
     * runtime exception.
     */
    private void RangeCheck(int index) {
    	if (Debug.ON && (index >= size || index < 0))
    		throw new IndexOutOfBoundsException(
    				"Index: "+index+", Size: "+size);
    }

    /**
     * Save the state of the ArrayList instance to a stream (that
     * is, serialize it).
     *
     * @serialData The length of the array backing the ArrayList
     *             instance is emitted (int), followed by all of its elements
     *             (each an Atom) in the proper order.
     */
    private synchronized void writeObject(java.io.ObjectOutputStream s)
        throws java.io.IOException{
    	// Write out element count, and any hidden stuff
    	s.defaultWriteObject();

        // Write out array length
        s.writeInt(elementData.length);

        // Write out all elements in the proper order.
        for (int i=0; i<size; i++)
            s.writeObject(elementData[i]);
    }

    /**
     * Reconstitute the ArrayList instance from a stream (that is,
     * deserialize it).
     */
    private synchronized void readObject(java.io.ObjectInputStream s)
        throws java.io.IOException, ClassNotFoundException {
    	// Read in size, and any hidden stuff
    	s.defaultReadObject();

        // Read in array length and allocate array
        int arrayLength = s.readInt();
        elementData = new Atom[arrayLength];

        // Read in all elements in the proper order.
        for (int i=0; i<size; i++)
            elementData[i] = (Atom)s.readObject();
    }
    
     /**
      * Returns an iterator over the elements in this list in proper
      * sequence. 
      *
      * This implementation returns a straightforward implementation of the
      * iterator interface, relying on the backing list's size(),
      * get(int), and remove(int) methods.
      *
      * Note that the iterator returned by this method will throw an
      * UnsupportedOperationException in response to its
      * remove method unless the list's remove(int) method is
      * overridden.
      *
      * This implementation can be made to throw runtime exceptions in the face
      * of concurrent modification, as described in the specification for the
      * (protected) modCount field.
      *
      * @return an iterator over the elements in this list in proper sequence.
      * 
      * @see #modCount
      */
     public AtomIterator iterator() {
     	return new AtomIteratorArrayListSimple(this);
     }
     
     public String toString() {
         StringBuffer buffer = new StringBuffer();
         for(int i=0; i<size; i++) {
             buffer.append(elementData[i].toString());
             buffer.append(" ");
         }
         return buffer.toString();
     }
 
}

