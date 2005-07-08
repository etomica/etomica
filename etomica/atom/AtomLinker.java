package etomica.atom;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.Atom;
import etomica.Debug;

/**
 * Class for constructing linked lists of Atoms.
 *
 * @author David Kofke
 */
public class AtomLinker implements java.io.Serializable {
    public final Atom atom;
    //TODO make these private and access with final methods, but test if performance suffers from change
    public transient AtomLinker next, previous;
    
    /**
     * Constructor throws exception if given atom is null.  Only
     * AtomLink.Tab instances can have a null atom field.
     */
    public AtomLinker(Atom a) {
        if(a == null && !(this instanceof Tab))
            throw new IllegalArgumentException("Error: cannot create AtomLinker with null atom");
        atom = a; 
        next = previous = this;
    }
    
    /**
     * Disconnects this linker from the one's previous and next to it, and
     * puts them in sequence, repairing the hole.
     */
    public void remove() {
        if (Debug.ON && previous.next != this || next.previous != this) {
            System.out.println("in AtomLinker.remove, this="+this+" previous.next="+previous.next+" next.previous="+next.previous);
            System.out.println("in AtomLinker.remove, atom="+atom+" previous.next.atom="+previous.next.atom+" next.previous.atom="+next.previous.atom);
            throw new RuntimeException("stop it");
        }
	    previous.next = next;
	    next.previous = previous;
	    if (Debug.ON) next = previous = null;
    }
    
    /**
     * Adds this linker to the list at the position before the given linker,
     * and after the linker that originally preceded the given one.
     */
    public void addBefore(AtomLinker newNext) {
        if(Debug.ON && this == newNext) {
            throw new IllegalArgumentException("Illegal attempt to link a linker to itself.");
        }
        next = newNext;
        previous = newNext.previous;
        previous.next = this;
        newNext.previous = this;
	}

    private void writeObject(java.io.ObjectOutputStream out)
    throws IOException
    {
        out.defaultWriteObject();
    }
    private void readObject(java.io.ObjectInputStream in)
    throws IOException, ClassNotFoundException
    {
        in.defaultReadObject();
        // linker is initially not "in the list"
        next = this;
        previous = this;
    }

	/**
	 * Creates a new tab that is flagged as not a header.
	 * @return Tab
	 */
	public static Tab newTab(AtomList list, int type) {
		return new Tab(list, type);
	}
	/**
	 * Creates a new tab that is flagged as being a header of a list.
	 * @return Tab
	 */
	public static Tab newHeader(AtomList list) {
		return new Tab(list, Tab.HEADER_TAB);
	}

    /**
     * Linker that does not reference an atom, but is used to hold or index
     * a position in the linked list.  This is the only type of linker that
     * can have a null atom field.  Tabs also hold references to each other in
     * the list, so it is possible to just directly from one tab to the next
     * (or previous) via the nextTab and previousTab fields.
     */
    public static class Tab extends AtomLinker {
        public Tab nextTab, previousTab;
        public AtomList list;
        private static int lastType = 0;
        public final int type;
        public final static int HEADER_TAB = Tab.requestTabType();
        public final static int ANY_TAB = 0xFFFFFFFF;
        /**
         * Private constructor.  Use AtomLinker methods newTab or newHeader, as
         * appropriate, to make new instance.
         * @param isHeader
         */
        private Tab(AtomList list, int type) {
            super(null);
            nextTab = previousTab = this;
            this.list = list;
            this.type = type;
        }
        
        public static int requestTabType() {
        	if (lastType == 0) {
        		lastType = 1;
        		return lastType;
        	}
        	// would java allow 32?
        	if (lastType == 1<<31) {
        		throw new RuntimeException("Too many tab types.  You can only have 31");
        	}
        	lastType = lastType << 1;
        	return lastType;
        }
        
        public final boolean isHeader() {return type == HEADER_TAB;}
        
        /**
         * Removes references to previous/next tabs while removing linker from list.
         */
        public void remove() {
        	if(isHeader()) throw new RuntimeException("Illegal attempt to remove header from list");
            super.remove();
	        previousTab.nextTab = nextTab;
	        nextTab.previousTab = previousTab;
	        nextTab = previousTab = this;
	        list = null;
        }
        
        /**
         * Adds references to previous/next tabs while adding the linker to the position
         * behind the given one.  Given linker must be in the list for this tab, otherwise
         * an IllegalArgumentException is thrown and the state of this tab is unchanged.
         */
        public void addBefore(AtomLinker newNext) {
	        nextTab = findNextTab(newNext);
	        this.list = nextTab.list;
            super.addBefore(newNext);
 	        previousTab = nextTab.previousTab;
	        previousTab.nextTab = this;
	        nextTab.previousTab = this;
	    }
       /**
        * Finds and returns the first Tab linker beginning at or after the given linker.
        */
        private AtomLinker.Tab findNextTab(AtomLinker e) {
            while(e.atom != null) e = e.next;
            return (AtomLinker.Tab)e;
        }
            
    }//end of Tab
        
}//end of AtomLinker
