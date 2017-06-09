/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IPotentialMolecular;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;

public class MayerFunctionMolecularThreeBody extends MayerFunctionThreeBody {

    protected final IPotentialMolecular p3;

    public MayerFunctionMolecularThreeBody(IPotentialMolecular p3) {
        this.p3 = p3;
    }

    protected double energy(IMoleculeList molecules, double[] r2) {
        return p3.energy(molecules);
    }

    public void setBox(Box box) {
        p3.setBox(box);
        super.setBox(box);
    }
}
