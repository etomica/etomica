/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.units.Dimension;
import etomica.units.Temperature;

/**
 * Integrator implements the algorithm used to move the atoms around and
 * generate new configurations in one or more boxs. All integrator techniques,
 * such as molecular dynamics or Monte Carlo, are implemented via subclasses of
 * this Integrator class. The Integrator's activities are managed via the
 * actions of the governing Controller.
 * 
 * @author David Kofke and Andrew Schultz
 */

public abstract class IntegratorBox extends Integrator {

    protected Box box;
    protected double temperature;
    protected boolean isothermal = false;
    protected DataSourceScalar meterPE;
    protected double currentPotentialEnergy;
    protected final PotentialMaster potentialMaster;

    public IntegratorBox(PotentialMaster potentialMaster, double temperature) {
        super();
        this.potentialMaster = potentialMaster;
        if (potentialMaster != null) {
            meterPE = new MeterPotentialEnergy(potentialMaster);
        }
        setTemperature(temperature);
    }

    /**
     * @return Returns the PotentialMaster.
     */
    public PotentialMaster getPotentialMaster() {
        return potentialMaster;
    }

    /**
     * Defines the actions taken by the integrator to reset itself, such as
     * required if a perturbation is applied to the simulated box (e.g.,
     * addition or deletion of a molecule). Also invoked when the
     * integrator is started or initialized. This also recalculates the 
     * potential energy.
     */
    public void reset() {
        super.reset();
        if (meterPE != null) {
            currentPotentialEnergy = meterPE.getDataAsScalar();
            if (currentPotentialEnergy == Double.POSITIVE_INFINITY) {
                System.err.println("overlap in configuration for "+box+" when resetting integrator");
                PotentialCalculationEnergySum.debug = true;
                meterPE.getDataAsScalar();
                PotentialCalculationEnergySum.debug = false;
                throw new ConfigurationOverlapException(box);
            }
        }
    }

    /**
     * Sets the temperature for this integrator
     */
    public void setTemperature(double t) {
        temperature = t;
    }

    /**
     * @return the integrator's temperature
     */
    public final double getTemperature() {
        return temperature;
    }

    /**
     * @return the dimenension of temperature (TEMPERATURE)
     */ 
    public final Dimension getTemperatureDimension() {
        return Temperature.DIMENSION;
    }

    /**
     * @return the potential energy of the box handled by this integrator
     */
    public double getPotentialEnergy() {
        return currentPotentialEnergy;
    }
    
	public void setIsothermal(boolean b) {
		isothermal = b;
	}

	public boolean isIsothermal() {
		return isothermal;
	}

	/**
	 * Performs activities needed to set up integrator to work on given box.
	 * 
	 * @return true if the box was successfully added to the integrator; false
	 *         otherwise
	 */
	public void setBox(Box p) {
	    box = p;
	    if(meterPE instanceof MeterPotentialEnergy){
	        ((MeterPotentialEnergy)meterPE).setBox(p);
	    }
	}

    public Box getBox() {
        return box;
    }

    public void setMeterPotentialEnergy(DataSourceScalar mpe){
        meterPE = mpe;
    }

    public DataSourceScalar getMeterPotentialEnergy() {
        return meterPE;
    }

}
