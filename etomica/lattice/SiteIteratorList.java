package etomica.lattice;

/**
 * Site iterator that traverses the elements of an SiteList.
 * IteratorDirective direction can set whether iterator goes up or down list.
 * Can be set or as a neighbor iterator, in which case the first indicated site
 * is skipped.
 *
 * @author David Kofke
 */
public final class SiteIteratorList implements SiteIterator {
    
	private SiteLinker next, header;
    private SiteList list;
        
	public SiteIteratorList() {
	    this(SiteList.NULL);
	}
	public SiteIteratorList(SiteIterator iterator) {
	    this(new SiteList(iterator));
	}
	public SiteIteratorList(SiteList list) {
	    setBasis(list);
	}
	
	public void reset(int index) {
	    if (index < 0 || index > list.size())
		throw new IndexOutOfBoundsException("Index: "+index+
						    ", Size: "+list.size());
	    if (index < list.size()/2) {
		    next = header.next;
		    for (int nextIndex=0; nextIndex<index; nextIndex++) next = next.next;
	    } else {
		    next = header;
		    for (int nextIndex=list.size(); nextIndex>index; nextIndex--) next = next.previous;
	    }
	}

	public boolean hasNext() {return next != header;}
	
    public void setBasis(SiteList list) {
        this.list = list;
        if(list == null) {
            header = null;
            next = null;
            return;
        }
        header = list.header;
        reset();
    }
    
    public void allSites(SiteAction action){
        for (SiteLinker e = header.next; e != header; e = e.next) action.actionPerformed(e.site);
    }
    
    public Site first() {
        if(header == null || header.next == null) return null;
        else return header.next.site;
    }
    
    /**
     * Sets iterator so that it is ready to go up its entire list of iterates.
     */
    public void reset() {
        reset(header.next);
    }
        
    /**
     * Resets in reference to the given site linker.  The
     * given linker (or its site) is the next iterate.
     * Does not check that the linker is an iterate of this iterator.
     */
    public void reset(SiteLinker linker) {
        if(linker == header || linker == null) {
            next = header;
            return;
        }
        else next = linker;
    }
    
    /**
     * Resets in reference to the given site.  Finds the site in the list and
     * calls reset(SiteLinker) in reference to its linker.  If site is not in list,
     * hasNext is set to false.
     */
    public void reset(Site site) {
        reset(list.linker(site));
    }

    /**
     * Returns true if the given site is in the list of iterates, false otherwise.
     */
	public boolean contains(Site site){
        return list.contains(site);
	}
	
	/**
	 * Returns the total number of iterates that can be returned by this iterator, for
	 * its current list basis.
	 */
	public int size(){ return list.size();}

	    
    public Site next() {
        return nextLinker().site;
    }
    
    public SiteLinker nextLinker() {
        SiteLinker nextLinker = next;
        next = next.next;
        return nextLinker;
    }

}//end of SiteIteratorList

