/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.freeenergy.npath;

import etomica.api.IBox;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMaster;
import etomica.units.Null;

/**
 * Meter whose sole purpose in life is return du/dw.
 *
 * Created by andrew on 5/8/17.
 */
public class MeterDUDW extends DataSourceScalar {

    protected final PotentialCalculationDUDW pc;
    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective id;
    protected final IBox box;

    public MeterDUDW(PotentialMaster potentialMatser, IBox box) {
        super("DUDW", Null.DIMENSION);
        this.potentialMaster = potentialMatser;
        this.box = box;
        id = new IteratorDirective();
        pc = new PotentialCalculationDUDW();
    }

    @Override
    public double getDataAsScalar() {
        pc.reset();
        potentialMaster.calculate(box, id, pc);
        return pc.getSum();
    }
}
