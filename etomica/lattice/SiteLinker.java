package etomica.lattice;

/**
 * Class for constructing linked lists of Sites.
 * Each Linker points to one site and another Linker, the next one in the list.
 * Although each site has built-in ability to link to one next and one previous site, these
 * Linkers are needed to construct other lists of sites, particularly for neighbor lists.
 *
 * @author David Kofke
 */
public final class SiteLinker implements java.io.Serializable {
    public final Site site;
    public SiteLinker next = null, previous = null;
    //Constructor
    public SiteLinker(Site a, SiteLinker next, SiteLinker previous) {
        site = a; 
        this.next = next;
        this.previous = previous;
        if(next != null) next.previous = this;
        if(previous != null) previous.next = this;
    }
}
