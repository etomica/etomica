/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Vector;

public interface IPotentialTorque extends PotentialSoft {

    Vector[][] gradientAndTorque(IAtomList atoms);
}
