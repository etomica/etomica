/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;

/**
 * Interface for a Mayer f-function, which takes 3+ atoms and
 * returns exp(-(u(atoms)-u(pairs))/kT) - 1
 *
 * @author Andrew Schultz
 */
public interface MayerFunctionNonAdditive {

    /**
     * returns exp(-beta*(U - Upair))
     * 
     * r2 is listed in the order
     *   (0,1),(0,2)...(0,n-1),(1,2),(1,3)...(1,n-1)...(n-2,n-1)
     */
	public double f(IMoleculeList molecules, double[] r2, double beta);

	/**
     * returns exp(-beta*(U - Upair))
     * 
     * This method allows an implementation to operate on distances and indices
     * instead of molecules.  The arrays can be longer than needed for
     * nMolecules (extra elements will be ignored).
	 */
	public double f(IMoleculeList molecules, int nMolecules, int[] moleculeIndices, double[] r2, double beta);

	public void setBox(Box box);
}
