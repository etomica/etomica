package etomica;

/**
 * Class for constructing linked lists of Atoms.
 *
 * @author David Kofke
 */
public class AtomLinker implements java.io.Serializable {
    public final Atom atom;
    public AtomLinker next, previous;
    
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
	    previous.next = next;
	    next.previous = previous;
	    next = previous = this;
    }
    
    /**
     * Adds this linker to the list at the position before the given linker,
     * and after the linker that originally preceded the given one.
     */
    public void addBefore(AtomLinker newNext) {
        next = newNext;
        previous = newNext.previous;
        previous.next = this;
        newNext.previous = this;
	}
	
	/**
	 * Moves the linker from its current position to the one before the given linker.
	 */
	public void moveBefore(AtomLinker newNext) {
	    remove();
	    addBefore(newNext);
	}

    //will add a reservoir of linkers, so this is used as constructor for now
    public static AtomLinker makeLinker(Atom a) {
        return new AtomLinker(a);
    }
    
	/**
	 * Creates a new tab that is flagged as not a header.
	 * @return Tab
	 */
	public static Tab newTab() {
		return new Tab(false);
	}
	/**
	 * Creates a new tab that is flagged as being a header of a list.
	 * @return Tab
	 */
	public static Tab newHeader() {
		return new Tab(true);
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
        private final boolean isHeader;
        /**
         * Private constructor.  Use AtomLinker methods newTab or newHeader, as
         * appropriate, to make new instance.
         * @param isHeader
         */
        private Tab(boolean isHeader) {
            super(null);
            nextTab = previousTab = this;
            this.isHeader = isHeader;
        }
        
        public final boolean isHeader() {return isHeader;}
        
        /**
         * Removes references to previous/next tabs while removing linker from list.
         */
        public void remove() {
            super.remove();
	        previousTab.nextTab = nextTab;
	        nextTab.previousTab = previousTab;
	        nextTab = previousTab = this;
        }
        
        /**
         * Adds references to previous/next tabs while adding the linker to the position
         * behind the given one.
         */
        public void addBefore(AtomLinker newNext) {
            super.addBefore(newNext);
	        nextTab = findNextTab(newNext);
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
