/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
