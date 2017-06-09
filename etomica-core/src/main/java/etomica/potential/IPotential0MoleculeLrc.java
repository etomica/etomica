/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IMolecule;
import etomica.api.IPotentialMolecular;

/**
 * Interface for a long-range correction potential.
 *
 * @author Andrew Schultz
 */
public interface IPotential0MoleculeLrc extends IPotentialMolecular {

    /**
     * Informs the potential of a target atom.  Null target atom indicates no
     * target atom.
     */
    public void setTargetMolecule(IMolecule targetAtom);

}
