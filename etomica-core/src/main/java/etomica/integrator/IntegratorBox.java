/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.box.Box;
import etomica.meta.annotations.IgnoreProperty;
import etomica.potential.PotentialMaster;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialMaster;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Dimensioned;
import etomica.units.dimensions.Temperature;

import java.util.Objects;

/**
 * Superclass of integrators that acts on a single box.
 *
 * @author David Kofke and Andrew Schultz
 */

public abstract class IntegratorBox extends Integrator {

    protected final PotentialMaster potentialMaster;
    protected final Box box;
    protected double temperature;
    protected boolean isothermal = false;
    protected DataSourceScalar meterPE;
    protected double currentPotentialEnergy;

    /**
     *
     * @param potentialMaster PotentialMaster instance used to compute energy etc.
     * @param temperature used by integration algorithm and/or to initialize velocities
     */
    public IntegratorBox(PotentialMaster potentialMaster, double temperature, Box box) {
        super();
        this.box = Objects.requireNonNull(box);
        this.potentialMaster = potentialMaster;
        if (potentialMaster != null) {
            meterPE = new MeterPotentialEnergy(potentialMaster, box);
        }
        setTemperature(temperature);
    }

    /**
     * @return the PotentialMaster instance used to compute energy etc.
     */
    @IgnoreProperty
    public PotentialMaster getPotentialMaster() {
        return potentialMaster;
    }

    /**
     * Performs superclass reset actions and recalculated currentPotentialEnergy
     *
     * @throws ConfigurationOverlapException if energy of current configuration is infinite
     */
    public void reset() {
        super.reset();
        if (meterPE != null) {
            currentPotentialEnergy = meterPE.getDataAsScalar();
            if (currentPotentialEnergy == Double.POSITIVE_INFINITY) {
                System.err.println("overlap in configuration for " + box + " when resetting integrator");
                PotentialCalculationEnergySum.debug = true;
                meterPE.getDataAsScalar();
                PotentialCalculationEnergySum.debug = false;
                throw new ConfigurationOverlapException(box);
            }
        }
    }

    /**
     * @return the integrator's temperature
     */
    @Dimensioned(dimension = Temperature.class)
    public final double getTemperature() {
        return temperature;
    }

    /**
     * @param t the new temperature
     */
    public void setTemperature(double t) {
        if(t < 0) {
            throw new IllegalArgumentException("Temperature cannot be negative");
        }
        temperature = t;
    }

    /**
     * @return the dimenension of temperature (TEMPERATURE)
     */
    public final Dimension getTemperatureDimension() {
        return Temperature.DIMENSION;
    }

    /**
     * @return the current value of the potential energy of the box handled by this integrator
     */
    public double getPotentialEnergy() {
        return currentPotentialEnergy;
    }

    /**
     *
     * @return true if the Integrator samples according to a specified temperature
     */
    public boolean isIsothermal() {
        return isothermal;
    }

    /**
     *
     * @param b specifies whether the Integrator should (if true) sample according to a specified temperature
     */
    public void setIsothermal(boolean b) {
        isothermal = b;
    }

    /**
     *
     * @return the Box instance on which this integrator acts
     */
    public Box getBox() {
        return box;
    }

    /**
     * @param box the Box to set
     */
    protected void setBox(Box box) {}

    /**
     *
     * @return the meter used to compute the potential energy
     */
    public DataSourceScalar getMeterPotentialEnergy() {
        return meterPE;
    }

    /**
     *
     * @param mpe the new meter used to compute the potential energy
     */
    public void setMeterPotentialEnergy(DataSourceScalar mpe) {
        meterPE = mpe;
    }

}
