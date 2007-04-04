package etomica.integrator;

import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.units.Dimension;
import etomica.units.Temperature;

/**
 * Integrator implements the algorithm used to move the atoms around and
 * generate new configurations in one or more phases. All integrator techniques,
 * such as molecular dynamics or Monte Carlo, are implemented via subclasses of
 * this Integrator class. The Integrator's activities are managed via the
 * actions of the governing Controller.
 * 
 * @author David Kofke and Andrew Schultz
 */

public abstract class IntegratorPhase extends Integrator {

    protected Phase phase;
    protected double temperature;
    protected boolean isothermal = false;
    protected MeterPotentialEnergy meterPE;
    protected double currentPotentialEnergy;

    public IntegratorPhase(PotentialMaster potentialMaster, double temperature) {
        super(potentialMaster);
        meterPE = new MeterPotentialEnergy(potentialMaster);
        setTemperature(temperature);
    }


    /**
     * Defines the actions taken by the integrator to reset itself, such as
     * required if a perturbation is applied to the simulated phase (e.g.,
     * addition or deletion of a molecule). Also invoked when the
     * integrator is started or initialized. This also recalculates the 
     * potential energy.
     */
    public void reset() throws ConfigurationOverlapException {
        meterPE.setPhase(phase);
        currentPotentialEnergy = meterPE.getDataAsScalar();
        if (currentPotentialEnergy == Double.POSITIVE_INFINITY) {
            System.err.println("overlap in configuration for "+phase+" when resetting integrator");
            throw new ConfigurationOverlapException(phase);
        }
    }

    /**
     * Perform any action necessary when neighbor lists are updated 
     */
    public void neighborsUpdated() {}
    
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
     * @return the potential energy of the phase handled by this integrator
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
	 * Performs activities needed to set up integrator to work on given phase.
	 * 
	 * @return true if the phase was successfully added to the integrator; false
	 *         otherwise
	 */
	public void setPhase(Phase p) {
	    phase = p;
	}

    public Phase getPhase() {
        return phase;
    }
}
