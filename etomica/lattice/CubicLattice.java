/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;


/**
 * Lattice or crystal with cubic symmetry, such that its spatial size can be
 * characterized by a single value, the lattice constant.
 */

/*
 * History
 * Created on Jan 5, 2005 by kofke
 */
public interface CubicLattice extends SpaceLattice {
    
    /**
     * Sets the size of the unit cell, the length of each edge.
     * @param latticeConstant the unit cell size
     */
    public void setLatticeConstant(double latticeConstant);
    
    /**
     * Accesses the size of the unit cell, the length of each edge of a single cell.
     */
    public double getLatticeConstant();
    
}
