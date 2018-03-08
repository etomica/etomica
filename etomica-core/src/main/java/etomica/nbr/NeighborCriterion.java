/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;
import etomica.atom.IAtom;
import etomica.box.Box;

/**
 * Atom filter used to specify whether two atoms are considered neighbors,
 * for the purpose of tabulating neighbor lists.  Neighbors are indicated 
 * if the AtomsetFilter.accept() method returns true.  Adds neighbor management
 * methods to set the box where the criterion is being applied, and to
 * indicate if neighbor list for atom needs updating.
 */
public interface NeighborCriterion {

    boolean accept(IAtom atom1, IAtom atom2);
    
    /**
     * Indicates whether the neighbor list for the given atom should
     * be updated, according to this criterion.  
     * @return true if the atom's list should be updated.
     */
    boolean needUpdate(IAtom atom);
	
    /**
     * Specifies the box where the criterion is being applied.  Sometimes
     * needed if the criterion depends on features of the box, such as the
     * boundary.
     */
    void setBox(Box box);
	
    /**
     * Indicates whether the atom has changed (e.g. moved) by an amount 
     * that might have caused its neighbor list to be invalid.  If this
     * method returns true, a neighbor list failure may have introduced
     * errors in the calculation.
     */
    boolean unsafe();
	
    /**
     * Indicates to criterion that given atom's neighbor list has just been 
     * updated, and that properties (e.g., record of atom's position) 
     * used by needUpdate and unsafe() methods should be reset.  
     */
    void reset(IAtom atom);
}
