/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.IAtom;
import etomica.space.Vector;

public interface IPotentialField {

    /**
     * Returns the energy between IAtom atom and the field.
     */
    double u(IAtom atom);

    /**
     * Computes the force (and adds it to f) for IAtom atom and returns the
     * energy due to the field.
     */
    double udu(IAtom atom, Vector f);

    default double uduTorque(IAtom atom, Vector f, Vector t) {
        return udu(atom, f);
    }

}
