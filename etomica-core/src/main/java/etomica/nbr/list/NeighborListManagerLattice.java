/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.math.function.IFunction;
import etomica.potential.BondingInfo;
import etomica.potential.IPotentialAtomic;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

public class NeighborListManagerLattice extends NeighborListManagerFasterer {
    public NeighborListManagerLattice(SpeciesManager sm, Box box, int cellRange, double nbrRange, BondingInfo bondingInfo) {
        super(sm, box, cellRange, nbrRange, bondingInfo);
    }

    protected int checkNbrPair(int i, int j, IAtom iAtom, IAtom jAtom, double rc2, Vector jbo, IPotentialAtomic[] iPotentials) {
        if (iPotentials[jAtom.getType().getIndex()] == null) return 0;

        Vector dr = space.makeVector();
        Vector ri = iAtom.getPosition();
        Vector rj = jAtom.getPosition();
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        dr.map(new IFunction() {
            @Override
            public double f(double x) {
                return Math.abs(x);
            }
        });
        int notZero = 0;
        for (int k=0; k<dr.getD(); k++) {
            if (dr.getX(k) > nbrRange) return 0;
            if (dr.getX(k) != 0) notZero++;
        }
        if (notZero > 1) return 0;
        return addAsNbrPair(i, j, iAtom, jAtom, jbo, iPotentials, dr);
    }

}
