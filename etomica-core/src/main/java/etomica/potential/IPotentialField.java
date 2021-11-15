/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;

public interface IPotentialField {

    /**
     * Returns the energy between IAtom atom and the field.
     */
    default double u(IAtom atom) {
        return energy(new AtomSetSinglet(atom));
    }

    /**
     * Computes the force (and adds it to f) for IAtom atom and returns the
     * energy due to the field.
     */
    double udu(IAtom atom, Vector f);

    default double uduTorque(IAtom atom, Vector f, Vector t) {
        return udu(atom, f);
    }

    default double energy(IAtomList atom) {
        return u(atom.get(0));
    }
}
