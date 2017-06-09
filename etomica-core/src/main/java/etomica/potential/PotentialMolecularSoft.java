/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.molecule.IMoleculeList;
import etomica.space.Tensor;
import etomica.space.Vector;


/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */
public interface PotentialMolecularSoft extends IPotentialMolecular {

    public double virial(IMoleculeList atoms);

    /**
	 * Returns the gradient of the potential as it applies to each atom in the 
     * given AtomSet, indicating how the energy would change as the position of 
     * the first atom is varied.  The method is allowed to return an array of
     * Vectors with fewer elements than the number of molecules in the
     * IMoleculeList.
	 */
	public Vector[] gradient(IMoleculeList molecules);
    
    /**
     * Returns the same gradient as gradient(IMoleculeList) and also adds in
     * the contribution of the molecules to the pressureTensor.  Their
     * contribution is added to the given Tensor.  This combined method exists
     * for computational efficiency.  Calculating the pressureTensor is
     * generally trivial once the gradient is known but often requires
     * intermediate information.
     */
    public Vector[] gradient(IMoleculeList molecules, Tensor pressureTensor);

}
