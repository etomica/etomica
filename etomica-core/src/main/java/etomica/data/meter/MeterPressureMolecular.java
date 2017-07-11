/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorBox;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationMolecularVirialSum;
import etomica.space.Space;
import etomica.units.dimensions.Pressure;

/**
 * Meter for evaluation of the soft-potential (molecular) pressure in a box.
 * Requires that temperature be set in order to calculation ideal-gas
 * contribution to pressure; default is to use zero temperature, which
 * causes this contribution to be omitted.
 *
 * @author Tai Boon Tan
 */
 
public class MeterPressureMolecular extends DataSourceScalar {
    
    public MeterPressureMolecular(Space space) {
    	super("Pressure",Pressure.dimension(space.D()));
    	dim = space.D();
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = true;
        virial = new PotentialCalculationMolecularVirialSum();
    }

    /**
     * Sets the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and integrator's 
     * temperature is used for the ideal gas contribution.
     */
    public void setIntegrator(IntegratorBox newIntegrator) {
        integrator = newIntegrator;
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
        if (integrator == null) {
            throw new IllegalStateException("You must call setIntegrator before using this class");
        }
    	virial.zeroSum();
        Box box = integrator.getBox();
        integrator.getPotentialMaster().calculate(box, iteratorDirective, virial);
        return (box.getMoleculeList().getMoleculeCount() / box.getBoundary().volume())*integrator.getTemperature() - virial.getSum()/(box.getBoundary().volume()*dim);
    }

    private static final long serialVersionUID = 1L;
    private IntegratorBox integrator;
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationMolecularVirialSum virial;
    private final int dim;
}
