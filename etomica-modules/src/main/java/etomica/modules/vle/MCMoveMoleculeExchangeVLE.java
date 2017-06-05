/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.vle;

import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.data.meter.MeterDensity;
import etomica.integrator.IntegratorBox;
import etomica.integrator.mcmove.MCMoveMoleculeExchange;
import etomica.space.Space;

public class MCMoveMoleculeExchangeVLE extends MCMoveMoleculeExchange {

    public MCMoveMoleculeExchangeVLE(PotentialMaster potentialMaster, IRandom random,
                                     Space space,
                                     IntegratorBox integrator1, IntegratorBox integrator2) {
        super(potentialMaster, random, space, integrator1, integrator2);
        meterDensity = new MeterDensity(space);
    }
    
    public boolean doTrial() {
        boolean success = super.doTrial();
        if (!success) return success;
        meterDensity.setBox(box1);
        double density1 = meterDensity.getDataAsScalar();
        meterDensity.setBox(box2);
        double density2 = meterDensity.getDataAsScalar();
        if (density2 > density1) {
            rejectNotify();
            return false;
        }
        return true;
    }
    
    private static final long serialVersionUID = 1L;
    protected double density;
    protected MeterDensity meterDensity;
}
