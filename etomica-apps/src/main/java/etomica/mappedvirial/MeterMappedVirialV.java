/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedvirial;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Pressure;

public class MeterMappedVirialV extends DataSourceScalar {

    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected final PotentialCalculationForceSum pcForce;
    protected final Box box;
    protected final IteratorDirective allAtoms;
    protected final AtomLeafAgentManager<Vector> forceManager;
    protected final PotentialCalculationMappedVirialV pc;
    
    public MeterMappedVirialV(Space space, PotentialMaster potentialMaster, Box box, int nbins) {
        super("pma",Pressure.DIMENSION);
        this.space = space;
        this.box = box;
        this.potentialMaster = potentialMaster;
        pcForce = new PotentialCalculationForceSum();
        if (box != null) {
            forceManager = new AtomLeafAgentManager<>(a -> space.makeVector(), box);
            pcForce.setAgentManager(forceManager);
        }
        else {
            forceManager = null;
        }
        pc = new PotentialCalculationMappedVirialV(space, box, nbins, forceManager);
        allAtoms = new IteratorDirective();
    }

    public PotentialCalculationMappedVirialV getPotentialCalculation() {
        return pc;
    }

    public double getDataAsScalar() {
        pcForce.reset();
        potentialMaster.calculate(box, allAtoms, pcForce);
        pc.reset();
        potentialMaster.calculate(box, allAtoms, pc);
        return pc.getPressure();
    }
}
