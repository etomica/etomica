/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IMoleculeList;
import etomica.space.Vector;

public interface IPotentialMolecularTorque extends PotentialMolecularSoft {

    public Vector[][] gradientAndTorque(IMoleculeList molecules);
}
