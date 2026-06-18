/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.amoeba;

import etomica.molecule.IMolecule;
import etomica.potential.IPotentialMolecular;

public interface IPotentialMoleculeSingle extends IPotentialMolecular {

    double energy(IMolecule molecule);

}
