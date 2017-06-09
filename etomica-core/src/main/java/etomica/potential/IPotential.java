/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;


import etomica.box.Box;

public interface IPotential {

    /**
     * Returns the range over which the potential applies.  IAtoms with a
     * greater separation do not interact.
     */
    public double getRange();

    /**
     * Informs the potential of the box on which it acts so that it can
     * properly consider the boundaries.
     */
    public void setBox(Box box);

    /**
     * The number of atoms on which the potential depends.
     */
    public int nBody();
}
