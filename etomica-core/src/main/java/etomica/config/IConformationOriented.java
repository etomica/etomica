/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.atom.IAtomList;
import etomica.space.IOrientation;

public interface IConformationOriented extends IConformation {

    /**
     * Assign coordinates to the atoms in the given list using the
     * given orientation for the molecule.
     */
    public void initializePositions(IAtomList atomList, IOrientation orientation);

}
