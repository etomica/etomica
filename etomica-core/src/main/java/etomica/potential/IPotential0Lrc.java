/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.api.IMolecule;
import etomica.api.IPotential;

/**
 * Interface for a long-range correction potential.
 *
 * @author Andrew Schultz
 */
public interface IPotential0Lrc extends IPotential {

    /**
     * Informs the potential of a target atom.  Null target atom indicates no
     * target atom.
     */
    public void setTargetAtom(IAtom targetAtom);

    /**
     * Informs the potential of a target atom.  Null target atom indicates no
     * target atom.
     */
    public void setTargetMolecule(IMolecule targetMolecule);
}
