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
    protected AtomLinker(Atom a) {
        if(a == null && !(this instanceof Tab))
            throw new IllegalArgumentException("Error: cannot create AtomLinker with null atom");
        atom = a; 
        next = previous = this;
    }
    
    public void remove() {
	    previous.next = next;
	    next.previous = previous;
	    next = previous = this;
    }
        
    public void addBefore(AtomLinker newNext) {
        next = newNext;
        previous = newNext.previous;
        previous.next = this;
        newNext.previous = this;
	}
	public void moveBefore(AtomLinker newNext) {
	    remove();
	    addBefore(newNext);
	}

    //will add a reservoir of linkers, so this is used as constructor for now
    public static AtomLinker makeLinker(Atom a) {
        return new AtomLinker(a);
    }
    
    public static class Tab extends AtomLinker {
        public Tab nextTab, previousTab;
        public Tab() {
            super(null);
            nextTab = previousTab = this;
        }
        
        public void remove() {
            super.remove();
	        previousTab.nextTab = nextTab;
	        nextTab.previousTab = previousTab;
	        nextTab = previousTab = this;
        }
        
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
            
    }
        
}//end of AtomLinker
