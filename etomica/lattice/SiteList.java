package etomica.lattice;
import etomica.Simulation;
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

public class SiteList implements java.io.Serializable
{
    public final SiteLinker header = new SiteLinker(null, null, null);
    private int size = 0;
    
    /**
     * Constructs an empty list.
     */
    public SiteList() {
        header.next = header.previous = header;
    }

    /**
     * Constructs a list containing the elements of the specified
     * collection, in the order they are returned by the collection's
     * iterator.
     */
     public SiteList(SiteIterator iterator) {
	    this();
	    addAll(iterator);
     }
     
     /**
      * Returns a randomly selected site from the list.
      */
     public Site getRandom() {
        if(size==0) return null;
        return entry((int)(Simulation.random.nextDouble()*size)).site;
     }

    /**
     * Returns the first element in this list.
     *
     * @return the first element in this list.
     */
    public Site getFirst() {
	    if (size==0) return null;
	    return header.next.site;
    }

    /**
     * Returns the last element in this list.
     *
     * @return the last element in this list.
     * @throws    NoSuchElementException if this list is empty.
     */
    public Site getLast()  {
	    if (size==0) return null;
	    return header.previous.site;
    }

    /**
     * Removes and returns the first element from this list.
     *
     * @return the first element from this list.
     * @throws    NoSuchElementException if this list is empty.
     */
    public Site removeFirst() {
        if(size == 0) return null;
	    Site first = header.next.site;
	    remove(header.next);
	    return first;
    }

    /**
     * Removes and returns the last element from this list.
     *
     * @return the last element from this list, or null if it is empty
     */
    public Site removeLast() {
        if(size == 0) return null;
	    Site last = header.previous.site;
	    remove(header.previous);
	    return last;
    }

    /**
     * Inserts the given element at the beginning of this list.
     */
    public void addFirst(Site o) {
	    addBefore(o, header.next);
    }

    /**
     * Appends the given element to the end of this list.  (Identical in
     * function to the <tt>add</tt> method; included only for consistency.)
     */
    public void addLast(Site o) {
	    addBefore(o, header);
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
    public boolean contains(Site o) {
        return indexOf(o) != -1;
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
     * @param o element to be appended to this list.
     * @return <tt>true</tt> (as per the general contract of
     * <tt>Collection.add</tt>).
     */
    public boolean add(Site o) {
	    addBefore(o, header);
        return true;
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
    public boolean remove(Site o) {
        if (o==null) return false;
        else {
            for (SiteLinker e = header.next; e != header; e = e.next) {
                if (o == e.site) {
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
    public boolean addAll(SiteIterator iterator) {
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
    public boolean addAll(int index, SiteIterator iterator) {
	    int numNew = iterator.size();
        if (numNew==0)
            return false;

        SiteLinker successor = (index==size ? header : entry(index));
        SiteLinker predecessor = successor.previous;
	    iterator.reset();
	    for (int i=0; i<numNew; i++) {
           // SiteLinker e = new SiteLinker.makeLinker(it.next(), successor, predecessor);
            SiteLinker e = new SiteLinker(iterator.next(), successor, predecessor);
            predecessor.next = e;
            predecessor = e;
        }
        successor.previous = predecessor;

        size += numNew;
        return true;
    }

    /**
     * Removes all of the elements from this list.
     */
    public void clear() {
        header.next = header.previous = header;
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
    public Site get(int index) {
        return entry(index).site;
    }

    /**
     * Replaces the element at the specified position in this list with the
     * specified element.
     *
     * @param index index of element to replace.
     * @param element element to be stored at the specified position.
     * @return the element previously at the specified position.
     * @throws IndexOutOfBoundsException if the specified index is out of
     *		  range (<tt>index &lt; 0 || index &gt;= size()</tt>).
     */
    public Site set(int index, Site element) {
        SiteLinker e = entry(index);
        new SiteLinker(element, e.next, e.previous);
        return e.site;
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
    public void add(int index, Site element) {
        addBefore(element, (index==size ? header : entry(index)));
    }

    /**
     * Removes the element at the specified position in this list.  Shifts any
     * subsequent elements to the left (subtracts one from their indices).
     * Returns the element that was removed from the list.
     *
     * @param index the index of the element to removed.
     * @return the element previously at the specified position.
     * 
     * @throws IndexOutOfBoundsException if the specified index is out of
     * 		  range (<tt>index &lt; 0 || index &gt;= size()</tt>).
     */
    public Site remove(int index) {
        SiteLinker e = entry(index);
        remove(e);
        return e.site;
    }

    /**
     * Return the indexed entry.
     */
    private SiteLinker entry(int index) {
        if (index < 0 || index >= size)
            throw new IndexOutOfBoundsException("Index: "+index+
                                                ", Size: "+size);
        SiteLinker e = header;
        if (index < size/2) {
            for (int i = 0; i <= index; i++)
                e = e.next;
        } else {
            for (int i = size; i > index; i--)
                e = e.previous;
        }
        return e;
    }


    // Search Operations
    
    /**
     * Returns the linker associated with the given site in this list.
     * Returns null if the site is not in the list.
     */
    public SiteLinker linker(Site site) {
        if(site == null) return null;
        for (SiteLinker e = header.next; e != header; e = e.next) {
            if (site.equals(e.site)) return e;
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
    public int indexOf(Site o) {
        if(o == null) return -1;
        int index = 0;
        for (SiteLinker e = header.next; e != header; e = e.next) {
            if (o.equals(e.site)) return index;
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
     * 	       specified element, or -1 if the list does not contain this
     * 	       element.
     */
    public int lastIndexOf(Site o) {
        if(o == null) return -1;
        int index = size;
        for (SiteLinker e = header.previous; e != header; e = e.previous) {
            index--;
            if (o.equals(e.site))
                return index;
        }
        return -1;
    }

    private SiteLinker addBefore(Site o, SiteLinker e) {
	    SiteLinker newSiteLinker = new SiteLinker(o, e, e.previous);
	    newSiteLinker.previous.next = newSiteLinker;
	    newSiteLinker.next.previous = newSiteLinker;
	    size++;
	    return newSiteLinker;
    }

    private void remove(SiteLinker e) {
	    if (e == header)
	        throw new NoSuchElementException();

	    e.previous.next = e.next;
	    e.next.previous = e.previous;
	    size--;
    }

    /**
     * Returns an array containing all of the elements in this list
     * in the correct order.
     *
     * @return an array containing all of the elements in this list
     * 	       in the correct order.
     */
    public Site[] toArray() {
	    Site[] result = new Site[size];
        int i = 0;
        for (SiteLinker e = header.next; e != header; e = e.next) result[i++] = e.site;
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
        SiteIterator iterator = phase.iteratorFactory().makeSiteIterator();
        iterator.reset();int i = 0;
        while(iterator.hasNext()) {
            Site next = iterator.next();
            next.parentGroup().setIndex(i);
            //System.out.println(next.parentGroup().index());
            i++;
        }
        SiteList[] list1 = new SiteList[8];
       // Site site = phase.lastSite().previousSite();System.out.println(site.toString());
        SiteIterator iterator1 = phase.iteratorFactory().makeSiteIterator();
        iterator1.reset();
        int ii = 0;
       while(iterator1.hasNext()&& ii < 8){
        list1[ii] = new SiteList();
        Site site = iterator1.next();System.out.println(site.parentGroup().index());
        list1[ii].addFirst(site);
        IteratorDirective directive = new IteratorDirective();
        iterator.reset(directive.set(site).set(IteratorDirective.BOTH));
        while(iterator.hasNext()) {
            Site next = iterator.next();
            if(next == site.nextSite()||next == site.previousSite()){
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
    public static final SiteList NULL = new SiteList();
    
}
