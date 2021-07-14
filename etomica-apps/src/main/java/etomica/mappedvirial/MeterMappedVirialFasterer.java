/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedvirial;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Pressure;

public class MeterMappedVirialFasterer extends DataSourceScalar {

    protected final PotentialCompute potentialMaster;
    protected final Box box;
    protected final PotentialCallbackMappedVirial pc;

    public MeterMappedVirialFasterer(PotentialCompute potentialMaster, Box box, int nbins) {
        super("pma",Pressure.DIMENSION);
        this.box = box;
        this.potentialMaster = potentialMaster;
        pc = new PotentialCallbackMappedVirial(box, potentialMaster, nbins);
    }

    public PotentialCallbackMappedVirial getPotentialCallback() {
        return pc;
    }

    public double getDataAsScalar() {
        potentialMaster.computeAll(true);
        pc.reset();
        potentialMaster.computeAll(false, pc);
        return pc.getPressure();
    }
}
