/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedvirial;

import etomica.atom.AtomLeafAgentManager;
import etomica.box.Box;
import etomica.box.storage.Tokens;
import etomica.box.storage.VectorStorage;
import etomica.data.DataSourceScalar;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Energy;

public class MeterMappedU extends DataSourceScalar {

    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected final PotentialCalculationForceSum pcForce;
    protected final Box box;
    protected final IteratorDirective allAtoms;
    protected final PotentialCalculationMappedEnergy pc;
    private final VectorStorage forces;

    public MeterMappedU(Space space, PotentialMaster potentialMaster, Box box, int nbins) {
        super("pma",Energy.DIMENSION);
        this.space = space;
        this.box = box;
        this.potentialMaster = potentialMaster;
        this.forces = box.getAtomStorage(Tokens.vectorsDefault());
        pcForce = new PotentialCalculationForceSum(forces);
        pc = new PotentialCalculationMappedEnergy(space, box, nbins, forces);
        allAtoms = new IteratorDirective();
    }

    public PotentialCalculationMappedEnergy getPotentialCalculation() {
        return pc;
    }

    public double getDataAsScalar() {
        pcForce.reset();
        potentialMaster.calculate(box, allAtoms, pcForce);
        pc.reset();
        potentialMaster.calculate(box, allAtoms, pc);
        return pc.getEnergy();
    }
}
