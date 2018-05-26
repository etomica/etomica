/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

/**
 * Interface for a lattice that has a finite, adjustable size.  Size is
 * quantified by the maximum value of the index in each dimension.  Provides
 * a method that permits all sites of the lattice to be returned in a 1-D array
 * of objects.
 */

/*
 * History
 * Created on Jan 3, 2005 by kofke
 */
public interface FiniteLattice<T> extends AbstractLattice<T> {

    /**
     * Returns an array containing all the sites in the lattice.
     */
    public T[] sites();

    /**
     * Returns an array with elements giving the maximum number of
     * index values in each dimension.  If size[i] = k, then permissible
     * values of index[i] are in the range (0, k-1) inclusive, where index is the 
     * index array given to the site method of AbstractLattice.
     */
    public int[] getSize();
    
    /**
     * Sets the array size, defined as described in the getSize method.
     */
    public void setSize(int[] size);
}
