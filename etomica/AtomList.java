package etomica;
import java.util.*;

/**
 * Linked list implementation of the <tt>List</tt> interface.  Implements all
 * optional list operations, and permits all elements (including
 * <tt>null</tt>).  In addition to implementing the <tt>List</tt> interface,
 * the <tt>LinkedList</tt> class provides uniformly named methods to
 * <tt>get</tt>, <tt>remove</tt> and <tt>insert</tt> an element at the
 * beginning and end of the list.  These operations allow linked lists to be
 * used as a stack, queue, or double-ended queue (deque).<p>
 *
 * All of the stack/queue/deque operations could be easily recast in terms of
 * the standard list operations.  They're included here primarily for
 * convenience, though they may run slightly faster than the equivalent List
 * operations.<p>
 *
 * All of the operations perform as could be expected for a doubly-linked
 * list.  Operations that index into the list will traverse the list from
 * the begining or the end, whichever is closer to the specified index.<p>
 *
 * <b>Note that this implementation is not synchronized.</b> If multiple
 * threads access a list concurrently, and at least one of the threads
 * modifies the list structurally, it <i>must</i> be synchronized
 * externally.  (A structural modification is any operation that adds or
 * deletes one or more elements; merely setting the value of an element is not
 * a structural modification.)  This is typically accomplished by
 * synchronizing on some object that naturally encapsulates the list.  If no
 * such object exists, the list should be "wrapped" using the
 * Collections.synchronizedList method.  This is best done at creation time,
 * to prevent accidental unsynchronized access to the list: <pre>
 *     List list = Collections.synchronizedList(new LinkedList(...));
 * </pre><p>
 *
 * The iterators returned by the this class's <tt>iterator</tt> and
 * <tt>listIterator</tt> methods are <i>fail-fast</i>: if the list is
 * structurally modified at any time after the iterator is created, in any way
 * except through the Iterator's own <tt>remove</tt> or <tt>add</tt> methods,
 * the iterator will throw a <tt>ConcurrentModificationException</tt>.  Thus,
 * in the face of concurrent modification, the iterator fails quickly and
 * cleanly, rather than risking arbitrary, non-deterministic behavior at an
 * undetermined time in the future.
 *
 * @author  Josh Bloch
 * @version 1.26 04/22/99
 * @see	    List
 * @see	    ArrayList
 * @see	    Vector
 * @see	    Collections#synchronizedList(List)
 * @since JDK1.2
 */

public class AtomList implements java.io.Serializable
{
    public final AtomLinker.Tab header = AtomLinker.newHeader();//modification for tab entry
    private int size = 0;
    
    /**
     * Constructs an empty list.
     */
    public AtomList() {
        header.next = header.previous = header;
        header.nextTab = header.previousTab = header;
    }

    /**
     * Constructs a list containing the elements of the specified
     * collection, in the order they are returned by the collection's
     * iterator.
     */
     public AtomList(AtomIterator iterator) {
	    this();
	    addAll(iterator);
     }
     
     /**
      * Makes a simple AtomLinker via the AtomLinker.makeLinker method.
      * May be overridden in subclasses to use linkers with other features.
      */
     protected AtomLinker makeLinker(Atom atom) {
        return AtomLinker.makeLinker(atom);
     }
     
	/**
	 * Positions the given tab in front of the header, effectively allowing it to serve
	 * as the header in this list.  It is assumed that the given tab lies in a 
	 * linked-list ring of atoms.  All atoms currently in this list are cleared and
	 * are automatically replaced by all atoms in new tab's ring.
	 * This is used for the cell neighbor list to place the neighborlisted atoms
	 * in a proper AtomList.  The integer argument newSize gives the number of
	 * (non-Tab) atoms in the pseudoHeader ring.
	 */
	 //used by IteratorFactoryCell.SequentialIterator
	public void setAsHeader(AtomLinker.Tab pseudoHeader, int newSize) {
	    header.remove();
	    header.addBefore(pseudoHeader);
	    size = newSize;
	//    size = 0;
	//    for(AtomLinker link=header.next; link!=header; link=link.next) {
	//        if(link.atom != null) size++;
	//    }
	}

     /**
      * Returns a randomly selected atom from the list.
      */
     public Atom getRandom() {
        if(size==0) return null;
        return entry((int)(Simulation.random.nextDouble()*size)).atom;
     }

    /**
     * Returns the first element in this list.
     *
     * @return the first element in this list.
     */
    public Atom getFirst() {
	    if (size==0) return null;
	    return firstEntry().atom;
    }

    /**
     * Returns the last element in this list.
     *
     * @return the last element in this list.
     * @throws    NoSuchElementException if this list is empty.
     */
    public Atom getLast()  {
	    if (size==0) return null;
	    return lastEntry().atom;
    }

    /**
     * Removes and returns the first (non-tab) element from this list.
     *
     * @return the first element from this list.
     * @throws    NoSuchElementException if this list is empty.
     */
    public Atom removeFirst() {
        if(size == 0) return null;
        AtomLinker firstEntry = firstEntry();
	    remove(firstEntry);
	    return firstEntry.atom;
    }

    /**
     * Removes and returns the last (non-tab) element from this list.
     *
     * @return the last element from this list, or null if it is empty
     */
    public Atom removeLast() {
        if(size == 0) return null;
	    AtomLinker lastEntry = lastEntry();
	    remove(lastEntry);
	    return lastEntry.atom;
    }

    /**
     * Inserts the given element at the beginning of this list.
     */
    public void addFirst(Atom atom) {
	    addBefore(atom, header.next);
    }
    
    /**
     * Inserts the given entry at the beginning of this list.
     */
    public void addFirst(AtomLinker linker) {
        addBefore(linker, header.next);
    }

    /**
     * Appends the given element to the end of this list.  (Identical in
     * function to the <tt>add</tt> method; included only for consistency.)
     */
    public void addLast(Atom atom) {
	    addBefore(atom, header);
    }

    /**
     * Returns <tt>true</tt> if this list contains the specified element.
     * More formally, returns <tt>true</tt> if and only if this list contains
     * at least one element <tt>e</tt> such that <tt>(o==null ? e==null
     * : o.equals(e))</tt>.
     *
     * @param o element whose presence in this list is to be tested.
     * @return <tt>true</tt> if this list contains the specified element.
     */
    public boolean contains(Atom atom) {
        return indexOf(atom) != -1;
    }

    /**
     * Returns the number of elements in this list.
     *
     * @return the number of elements in this list.
     */
    public int size() {
	    return size;
    }

    /**
     * Appends the specified element to the end of this list.
     *
     * @param atom atom to be appended to this list.
     * @return <tt>true</tt> (as per the general contract of
     * <tt>Collection.add</tt>).
     */
    public boolean add(Atom atom) {
	    addBefore(atom, header);
        return true;
    }
    
    /**
     * Appends the specified atom linker to the end of this list.
     */
    public void add(AtomLinker link) {
        addBefore(link, header);
    }

    /**
     * Removes the first occurrence of the specified element in this list.  If
     * the list does not contain the element, it is unchanged.  More formally,
     * removes the element with the lowest index <tt>i</tt> such that
     * <tt>(o==null ? get(i)==null : o.equals(get(i)))</tt> (if such an
     * element exists).
     *
     * @param o element to be removed from this list, if present.
     * @return <tt>true</tt> if the list contained the specified element.
     */
    public boolean remove(Atom o) {
        if (o==null) return false;
        else {
            for (AtomLinker e = header.next; e != header; e = e.next) {
                if (o == e.atom) {
                    remove(e);
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Appends all of the elements in the specified collection to the end of
     * this list, in the order that they are returned by the specified
     * collection's iterator.  The behavior of this operation is undefined if
     * the specified collection is modified while the operation is in
     * progress.  (This implies that the behavior of this call is undefined if
     * the specified Collection is this list, and this list is nonempty.)
     *
     * @param index index at which to insert first element
     *			  from the specified collection.
     * @param c elements to be inserted into this list.
     * 
     * @throws IndexOutOfBoundsException if the specified index is out of
     *         range (<tt>index &lt; 0 || index &gt; size()</tt>).
     */
    public boolean addAll(AtomIterator iterator) {
        return addAll(size, iterator);
    }

    /**
     * Inserts all of the elements in the specified collection into this
     * list, starting at the specified position.  Shifts the element
     * currently at that position (if any) and any subsequent elements to
     * the right (increases their indices).  The new elements will appear
     * in the list in the order that they are returned by the
     * specified collection's iterator.
     *
     * @param index index at which to insert first element
     *		    from the specified collection.
     * @param c elements to be inserted into this list.
     * @throws IndexOutOfBoundsException if the specified index is out of
     *            range (<tt>index &lt; 0 || index &gt; size()</tt>).
     */
    public boolean addAll(int index, AtomIterator iterator) {
	    int numNew = iterator.size();
        if (numNew==0)
            return false;

        AtomLinker successor = (index==size ? header : entry(index));
        AtomLinker predecessor = successor.previous;
	    iterator.reset();//should remove this
	    for (int i=0; i<numNew; i++) {
	        makeLinker(iterator.next()).addBefore(successor);
	    }

        size += numNew;
        return true;
    }

    /**
     * Removes all of the elements (including tabs) from this list.
     */
    public void clear() {
        header.remove(); //this has the effect of disconnecting the header from
                         //all other elements, thereby removing them from this list
	    size = 0;
    }


    // Positional Access Operations

    /**
     * Returns the element at the specified position in this list.
     *
     * @param index index of element to return.
     * @return the element at the specified position in this list.
     * 
     * @throws IndexOutOfBoundsException if the specified index is is out of
     * range (<tt>index &lt; 0 || index &gt;= size()</tt>).
     */
    public Atom get(int index) {
        return entry(index).atom;
    }

    /**
     * Replaces the atom at the specified position in this list with the
     * specified atom.
     *
     * @param index index of element to replace.
     * @param element atom to be stored at the specified position.
     * @return the atom previously at the specified position.
     * @throws IndexOutOfBoundsException if the specified index is out of
     *		  range (<tt>index &lt; 0 || index &gt;= size()</tt>).
     */
    public Atom set(int index, Atom atom) {
        AtomLinker e = entry(index);
        makeLinker(atom).addBefore(e);
        e.remove();
        return e.atom;
    }

    /**
     * Inserts the specified element at the specified position in this list.
     * Shifts the element currently at that position (if any) and any
     * subsequent elements to the right (adds one to their indices).
     *
     * @param index index at which the specified element is to be inserted.
     * @param element element to be inserted.
     * 
     * @throws IndexOutOfBoundsException if the specified index is out of
     *		  range (<tt>index &lt; 0 || index &gt; size()</tt>).
     */
    public void add(int index, Atom element) {
        addBefore(element, (index==size ? header : entry(index)));
    }

    /**
     * Removes the element at the specified position in this list.  
     * Returns the element that was removed from the list.
     *
     * @param index the index of the element to removed.
     * @return the element previously at the specified position.
     * 
     * @throws IndexOutOfBoundsException if the specified index is out of
     * 		  range (<tt>index &lt; 0 || index &gt;= size()</tt>).
     */
    public Atom remove(int index) {
        AtomLinker e = entry(index);
        remove(e);
        return e.atom;
    }

    /**
     * Return the indexed entry.
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
     *
     * @return the first entry in this list.
     */
    public AtomLinker firstEntry() {
	    if (size==0) return null;
	    AtomLinker entry = header.next;
	    while(entry.atom == null) entry = entry.next;//modification for tab entry
	    return entry;
    }

    /**
     * Returns the last entry (linker with non-null atom) in this list.
     *
     * @return the last entry in this list.
     * @throws    NoSuchElementException if this list is empty.
     */
    public AtomLinker lastEntry()  {
	    if (size==0) return null;
	    AtomLinker entry = header.previous;
	    while(entry.atom == null) entry = entry.previous;
	    return entry;
    }

    // Search Operations
    
    /**
     * Returns the linker associated with the given atom in this list.
     * Returns null if the atom is not in the list.
     */
    public AtomLinker entry(Atom atom) {
        if(atom == null) return null;
        for (AtomLinker e = header.next; e != header; e = e.next) {
            if (atom.equals(e.atom)) return e;
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
     * @param o element to search for.
     * @return the index in this list of the first occurrence of the
     * 	       specified element, or -1 if the list does not contain this
     * 	       element.
     */
    public int indexOf(Atom o) {
        if(o == null) return -1;
        int index = 0;
        for (AtomLinker e = header.next; e != header; e = e.next) {
            if (o.equals(e.atom)) return index;
            if(e.atom != null) index++;//modification for tab entry
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
     * 	       specified element, or -1 if the list does not contain this
     * 	       element.
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

    public AtomLinker addBefore(Atom atom, AtomLinker e) {
	    return addBefore(makeLinker(atom), e);
	}
	
	/**
	 * Places the linker given by the first argument before the linker
	 * given by the second one.   The first linker should not already
	 * be in the list, and the second one should, but no checks are made
	 * to ensure this.  All methods that result in the addition of one or
	 * more atoms/linkers to the list work through this method.
	 */
	public AtomLinker addBefore(AtomLinker newAtomLinker, AtomLinker e) {
	    if(newAtomLinker.atom != null) size++;//modification for tab entry
	    newAtomLinker.addBefore(e);	        
	    return newAtomLinker;
    }
    
    
    /**
     * Moves the first linker to the position before the second one.
     * Equivalent to remove(moving); addBefore(moving, newNext); but
     * saves incrementing and decrementing size.
     * Assumes both entries are already in the list, and that they are 
     * not the same entry.
     */
    public void moveBefore(AtomLinker moving, AtomLinker newNext) {
        moving.moveBefore(newNext);
    }        

    /**
     * Removes the given linker from this list.  Does not check that linker
     * is in fact contained in list.  All methods that remove an atom/linker
     * from this list work through this method.
     */
    public void remove(AtomLinker e) {
	    if (e.atom != null) size--;//decrement size counter if not removing a Tab
	    else if(e == header) throw new NoSuchElementException();
	    e.remove();
    }

    /**
     * Returns an array containing all of the elements in this list
     * in the correct order (and excluding all tabs).
     *
     * @return an array containing all of the elements in this list
     * 	       in the correct order.
     */
    public Atom[] toArray() {
	    Atom[] result = new Atom[size];
        int i = 0;
        for (AtomLinker e = header.next; e != header; e = e.next) {
            if(e.atom != null) result[i++] = e.atom;//modification for tab entry
        }
	    return result;
    }

    
    //main method to demonstrate and test this class
/*    public static void main(String[] args) throws java.io.IOException {
        
        java.io.BufferedReader in = new java.io.BufferedReader(new java.io.InputStreamReader(System.in));
        Simulation.instance = new Simulation(new Space1D());
        Phase phase = new Phase();
       // Phase phase2 = new Phase();
        Species species = new SpeciesSpheres();
        species.setNMolecules(8);
        Simulation.instance.elementCoordinator.go();
        System.out.println("Hello");
        AtomIterator iterator = phase.iteratorFactory().makeAtomIterator();
        iterator.reset();int i = 0;
        while(iterator.hasNext()) {
            Atom next = iterator.next();
            next.parentGroup().setIndex(i);
            //System.out.println(next.parentGroup().index());
            i++;
        }
        AtomList[] list1 = new AtomList[8];
       // Atom atom = phase.lastAtom().previousAtom();System.out.println(atom.toString());
        AtomIterator iterator1 = phase.iteratorFactory().makeAtomIterator();
        iterator1.reset();
        int ii = 0;
       while(iterator1.hasNext()&& ii < 8){
        list1[ii] = new AtomList();
        Atom atom = iterator1.next();System.out.println(atom.parentGroup().index());
        list1[ii].addFirst(atom);
        IteratorDirective directive = new IteratorDirective();
        iterator.reset(directive.set(atom).set(IteratorDirective.BOTH));
        while(iterator.hasNext()) {
            Atom next = iterator.next();
            if(next == atom.nextAtom()||next == atom.previousAtom()){
                list1[ii].addLast(next);System.out.println(next.parentGroup().index());
            }
        }
        if(list1[ii].size() == 2){
         System.out.println(ii+"\t"+list1[ii].size()+"\t"+list1[ii].get(0).toString()+"\t"+list1[ii].get(1).toString()+"\t");
        }else{
         System.out.println(ii+"\t"+list1[ii].size()+"\t"+list1[ii].get(0).toString()+"\t"+list1[ii].get(1).toString()+"\t"+list1[ii].get(2).toString());
        }
        i++;
       }
    }//end of main
   */ 
    //meant to be a permanently empty list that can be used as a place holder when needed
    public static final AtomList NULL = new AtomList();
    
}
