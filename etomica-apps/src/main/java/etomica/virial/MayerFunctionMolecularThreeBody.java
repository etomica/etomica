/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;

public class MayerFunctionMolecularThreeBody extends MayerFunctionThreeBody {

    protected final IPotentialMolecular p3;

    public MayerFunctionMolecularThreeBody(IPotentialMolecular p3) {
        this.p3 = p3;
    }

    protected double energy(IMoleculeList molecules, double[] r2) {
        return p3.energy(molecules);
    }

    protected double energy(IMoleculeList molecules, double rAB2, double rAC2, double rBC2) { return p3.energy(molecules);}

}
