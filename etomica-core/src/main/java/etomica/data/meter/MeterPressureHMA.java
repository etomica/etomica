/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorBox;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationVirialSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.units.dimensions.Pressure;

public class MeterPressureHMA extends DataSourceScalar {

    // Constructor takes a Space object.
    public MeterPressureHMA(Space space) {
        super("Pressure",Pressure.dimension(space.D()));    // passes 'label' and dimensions for pressure.
        dim = space.D();    // creates variable dim, which holds the dimensions of the space.
        iteratorDirective = new IteratorDirective();    // Contains instructions on how the atoms are selected.
        iteratorDirective.includeLrc = false;    // Indicates if long-range contributions should be included.
        virial = new PotentialCalculationVirialSum();       // Creates a PotentialCalculation object which will compute viral across all iterated atoms.

    }

    /**
     * Sets the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and integrator's 
     * temperature is used for the ideal gas contribution.  Alternatively, you
     * can set the potentialMaster, temperature and box separately.
     */
    public void setIntegrator(IntegratorBox newIntegrator) {
        integrator = newIntegrator;
    }       // Sets the integrator. e.g IntegratorMC()

    public void setPotentialMaster(PotentialMaster newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }   // This controls all potentials within a simulation. e.g. If you need to calculate the energy, forces, etc.

    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }

    public void setBox(Box newBox) {
        box = newBox;
    }

    /**
     * Returns the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and integrator's 
     * temperature is used for the ideal gas contribution.
     */
    public IntegratorBox getIntegrator() {
        return integrator;
    }

    /**
     * Sets flag indicating whether calculated energy should include
     * long-range correction for potential truncation (true) or not (false).
     */
    public void setIncludeLrc(boolean b) {
        iteratorDirective.includeLrc = b;
    }

    /**
     * Indicates whether calculated energy should include
     * long-range correction for potential truncation (true) or not (false).
     */
    public boolean isIncludeLrc() {
        return iteratorDirective.includeLrc;
    }

    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    public double getDataAsScalar() {
        if (integrator == null && (potentialMaster == null || box == null)) {
            throw new IllegalStateException("You must call setIntegrator before using this class");
        }
        virial.zeroSum();
        Box b = box;
        if (b == null) {
            b = integrator.getBox();
        }
        if (potentialMaster != null) {
            potentialMaster.calculate(b, iteratorDirective, virial);
        }
        else {
            integrator.getPotentialMaster().calculate(b, iteratorDirective, virial);
        }
        double temp = (integrator != null) ? integrator.getTemperature() : temperature;
        //System.out.println("fac="+(1/(box.getBoundary().volume()*box.getSpace().D())));
        return (b.getMoleculeList().size() / b.getBoundary().volume())*temp - virial.getSum()/(b.getBoundary().volume()*dim);
    }

    private IntegratorBox integrator;
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationVirialSum virial;
    protected PotentialMaster potentialMaster;
    protected double temperature;
    protected Box box;
    private final int dim;
}