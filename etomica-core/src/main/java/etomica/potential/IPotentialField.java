/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Tensor;
import etomica.space.Vector;

public interface IPotentialField extends PotentialSoft {

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
    default double udu(IAtom atom, Vector f) {
        AtomSetSinglet singlet = new AtomSetSinglet(atom);
        double u = energy(singlet);
        Vector[] g1 = gradient(singlet);
        f.ME(g1[0]);
        return u;
    }

    default int nBody() {
        return 1;
    }

    default double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    default void setBox(Box box) {
    }

    default double virial(IAtomList atom) {
        return 0;
    }

    default double energy(IAtomList atom) {
        return u(atom.get(0));
    }

    default Vector[] gradient(IAtomList atoms) {
        IAtom atom = atoms.get(0);
        Vector f = Vector.d(atom.getPosition().getD());
        udu(atom, f);
        f.TE(-1);
        return new Vector[]{f};
    }

    default Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }
}
