/*
 * History
 * Created on Dec 17, 2004 by kofke
 */
package etomica.lattice;


/**
 * Interface for a class that generates the sites populating a lattice.
 */
public interface SiteFactory {

    /**
     * Constructs a site for placement in a lattice.  Arguments
     * are not necessarily used to construct the site, but are provided
     * by the interface in case this information is to be held or used
     * by the site. 
     * @param lattice the lattice where the site will reside
     * @param coord the index by which the site can be accessed from the lattice. New site should
     * make its own copy of this if it uses the coord.
     */
    public Object makeSite(AbstractLattice lattice, int[] coord);
}
