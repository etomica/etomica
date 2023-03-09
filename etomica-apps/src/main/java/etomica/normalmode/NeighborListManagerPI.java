/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.nbr.list.NeighborListManager;
import etomica.potential.BondingInfo;
import etomica.potential.IPotential2;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

public class NeighborListManagerPI extends NeighborListManager {


    public NeighborListManagerPI(SpeciesManager sm, Box box, int cellRange, double nbrRange, BondingInfo bondingInfo) {
        super(sm, box, cellRange, nbrRange, bondingInfo);
    }

    protected int checkNbrPair(int i, int j, IAtom iAtom, IAtom jAtom, double rc2, Vector jbo, IPotential2[] iPotentials) {
        if (iAtom.getIndex() != jAtom.getIndex()) return 0;
        return super.checkNbrPair(i, j, iAtom, jAtom, rc2, jbo, iPotentials);
    }

}
