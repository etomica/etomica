/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.molecule;

import etomica.atom.ChargeSource;
import etomica.atom.DipoleSourceAtomic;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.space.Boundary;
import etomica.space.Vector;

public class DipoleSourceMolecularGeneric implements DipoleSourceMolecular {

    protected final Boundary boundary;
    protected final DipoleSourceAtomic dipoleSourceAtomic;
    protected final ChargeSource chargeSource;
    protected final Vector dipole;
    protected final Vector dr;

    public DipoleSourceMolecularGeneric(Box box, ChargeSource chargeSource, DipoleSourceAtomic dipoleSourceAtomic) {
        this.boundary = box.getBoundary();
        dipole = box.getSpace().makeVector();
        dr = box.getSpace().makeVector();
        this.chargeSource = chargeSource;
        this.dipoleSourceAtomic = dipoleSourceAtomic;
    }

    @Override
    public Vector getDipole(IMolecule molecule) {
        dipole.E(0);
        Vector com = MoleculePositionCOMPBC.com(boundary, molecule);
        for (IAtom a : molecule.getChildList()) {
            if (dipoleSourceAtomic != null) {
                dipole.PE(dipoleSourceAtomic.getDipole(a));
            }
            if (chargeSource != null) {
                dr.Ev1Mv2(a.getPosition(), com);
                boundary.nearestImage(dr);
                dipole.PEa1Tv1(chargeSource.getCharge(a.getType()), dr);
            }
        }
        return dipole;
    }
}
