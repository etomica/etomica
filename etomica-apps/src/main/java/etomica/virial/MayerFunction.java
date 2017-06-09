/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IPotential;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;

/**
 * @author kofke
 *
 * Abstract class for a Mayer f-function, which takes a pair of atoms and
 * returns exp(-u(pair)/kT) - 1
 */
public interface MayerFunction {

    /**
     * returns Mayer function between atoms in the pair at temperature
     * 1/beta
     * @param r2 TODO
     */
	public double f(IMoleculeList pair, double r2, double beta);

	/**
	 * @return
	 */
	public IPotential getPotential();
	
	public void setBox(Box box);
}
