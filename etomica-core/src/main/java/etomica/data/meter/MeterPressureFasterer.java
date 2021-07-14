/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.IDataSourcePotential;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Pressure;

/**
 * Meter for evaluation of the soft-potential pressure in a box.
 * Requires that temperature be set in order to calculation ideal-gas
 * contribution to pressure; default is to use zero temperature, which
 * causes this contribution to be omitted.
 *
 * @author Andrew Schultz
 */

public class MeterPressureFasterer extends DataSourceScalar implements IDataSourcePotential {

    protected final PotentialCompute potentialCompute;
    protected double temperature;
    protected final Box box;
    private final int dim;
    protected boolean callComputeAll = true;

    public MeterPressureFasterer(Box box, PotentialCompute potentialCompute) {
        super("Pressure", Pressure.dimension(box.getSpace().D()));
        dim = box.getSpace().D();
        this.box = box;
        this.potentialCompute = potentialCompute;
    }

    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }

    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    public double getDataAsScalar() {
        if (callComputeAll) potentialCompute.computeAll(potentialCompute.needForcesForVirial());
        //System.out.println("fac="+(1/(box.getBoundary().volume()*box.getSpace().D())));
        double vol = box.getBoundary().volume();
        return (box.getMoleculeList().size() / vol) * temperature
                - potentialCompute.getLastVirial() / (vol * dim);
    }

    @Override
    public void doCallComputeAll(boolean callComputeAll) {
        this.callComputeAll = callComputeAll;
    }
}
