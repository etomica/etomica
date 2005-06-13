package etomica.atom;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.Debug;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 12, 2005 by kofke
 */
public class AtomListTabbed extends AtomList {

    /**
     * 
     */
    public AtomListTabbed() {
        super();
    }

    /**
     * @param atoms
     */
    public AtomListTabbed(Atom[] atoms) {
        super(atoms);
    }

    /**
     * @param list
     */
    public AtomListTabbed(AtomList list) {
        super(list);
    }

    /**
     * @param iterator
     */
    public AtomListTabbed(AtomIterator iterator) {
        super(iterator);
    }

    /**
     * Removes all of the elements (including tabs) from this list.
     */
    public void clear() {
        header.nextTab = header.previousTab = header;
        header.previous = header.next = header;
        size = 0;
    }

    /**
     * Return the indexed entry, counting from 0.
     */
    public AtomLinker entry(int index) {
        if (index < 0 || index >= size)
            throw new IndexOutOfBoundsException("Index: "+index+
                                                ", Size: "+size);
        AtomLinker e = header;
        if (index < size/2) {
            for (int i = 0; i <= index; i++) {
                e = e.next;
                if(e.atom == null) i--;//modification for tab entry
            }
        } else {
            for (int i = size; i > index; i--) {
                e = e.previous;
                if(e.atom == null) i++;//modification for tab entry
            }
        }
        return e;
    }
    
    /**
     * Returns the first entry (linker with non-null atom) in this list.
     * Returns null if the list is empty.
     *
     * @return the first entry in this list.
     */
    public AtomLinker firstEntry() {
        if (isEmpty()) return null;
        AtomLinker entry = header.next;
        while(entry.atom == null) entry = entry.next;//modification for tab entry
        return entry;
    }

    /**
     * Returns the last entry (linker with non-null atom) in this list.
     * Returns null if the list is empty.
     */
    public AtomLinker lastEntry()  {
        if (isEmpty()) return null;
        AtomLinker entry = header.previous;
        while(entry.atom == null) entry = entry.previous;
        return entry;
    }

    /**
     * Returns the linker associated with the given atom in this list.
     * Returns null if the atom is not in the list, or if the atom is null.
     */
    public AtomLinker entry(Atom atom) {
        if(atom == null) return null;
        for (AtomLinker e = header.next; e != header; e = e.next) {
            if (e.atom!=null && atom.equals(e.atom)) return e;
        }
        return null;
    }

    /**
     * Returns the index in this list of the first occurrence of the
     * specified element, or -1 if the List does not contain this
     * element.  More formally, returns the lowest index i such that
     * <tt>(o==null ? get(i)==null : o.equals(get(i)))</tt>, or -1 if
     * there is no such index.
     *
     * @param atom element to search for.
     * @return the index in this list of the first occurrence of the
     *         specified element, or -1 if the list does not contain this
     *         element.
     */
    public int indexOf(Atom atom) {
        if(atom == null) return -1;
        int index = 0;
        for (AtomLinker e = header.next; e != header; e = e.next) {
            if (e.atom == null) continue; //skip tabs
            if (atom.equals(e.atom)) return index;
            index++;
        }
        return -1;
    }

    /**
     * Returns the index in this list of the last occurrence of the
     * specified element, or -1 if the list does not contain this
     * element.  More formally, returns the highest index i such that
     * <tt>(o==null ? get(i)==null : o.equals(get(i)))</tt>, or -1 if
     * there is no such index.
     *
     * @param o element to search for.
     * @return the index in this list of the last occurrence of the
     *         specified element, or -1 if the list does not contain this
     *         element.
     */
    public int lastIndexOf(Atom atom) {
        if(atom == null) return -1;
        int index = size;
        for (AtomLinker e = header.previous; e != header; e = e.previous) {
            if(e.atom != null) index--;//modification for tab entry
            if (atom.equals(e.atom))
                return index;
        }
        return -1;
    }

    /**
     * Places the linker given by the first argument before the linker
     * given by the second one.  Overrides parent class to permit addition
     * of Tab.
     */
    public AtomLinker addBefore(AtomLinker newAtomLinker, AtomLinker e) {
        if(newAtomLinker.atom != null) size++;//modification for tab entry
        newAtomLinker.addBefore(e);         
        return newAtomLinker;
    }

    /**
     * Removes the given linker from this list.  Does not check that linker
     * is in fact contained in list.  All methods that remove an atom/linker
     * from this list work through this method.
     */
    public void remove(AtomLinker e) {
        if(Debug.ON && !inList(e)) throw new IllegalArgumentException("Illegal attempt to remove a linker that is not in the list");
        if(e.atom != null) size--;//decrement size counter
        e.remove();//e.remove
    }

    /**
     * Returns an array containing all of the elements in this list
     * in the correct order (and excluding all tabs).
     *
     * @return an array containing all of the elements in this list
     *         in the correct order.
     */
    public Atom[] toArray() {
        Atom[] result = new Atom[size];
        int i = 0;
        for (AtomLinker e = header.next; e != header; e = e.next) {
            if(e.atom != null) result[i++] = e.atom;//modification for tab entry
        }
        return result;
    }

}
