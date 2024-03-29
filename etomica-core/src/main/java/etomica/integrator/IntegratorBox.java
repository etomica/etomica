/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.box.Box;
import etomica.exception.ConfigurationOverlapException;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Dimensioned;
import etomica.units.dimensions.Temperature;
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Objects;

/**
 * Superclass of integrators that acts on a single box.
 *
 * @author David Kofke and Andrew Schultz
 */

public abstract class IntegratorBox extends Integrator implements Statefull {

    protected final PotentialCompute potentialCompute;
    protected final Box box;
    protected final Space space;
    protected double temperature;
    protected boolean isothermal = false;
    protected double currentPotentialEnergy;

    /**
     * @param potentialCompute PotentialMaster instance used to compute energy etc.
     * @param temperature     used by integration algorithm and/or to initialize velocities
     */
    public IntegratorBox(PotentialCompute potentialCompute, double temperature, Box box) {
        super();
        this.box = Objects.requireNonNull(box);
        this.space = box.getSpace();
        this.potentialCompute = potentialCompute;
        if (potentialCompute != null) this.getEventManager().addListener(this.potentialCompute.makeIntegratorListener());
        setTemperature(temperature);
    }

    /**
     * @return the PotentialMaster instance used to compute energy etc.
     */
    public PotentialCompute getPotentialCompute() {
        return potentialCompute;
    }

    /**
     * Performs superclass reset actions and recalculated currentPotentialEnergy
     *
     * @throws ConfigurationOverlapException if energy of current configuration is infinite
     */
    public void reset() {
        super.reset();
        currentPotentialEnergy = potentialCompute.computeAll(false);
        if (currentPotentialEnergy == Double.POSITIVE_INFINITY) {
            System.err.println("overlap in configuration for " + box + " when resetting integrator");
            potentialCompute.computeAll(false);
            throw new ConfigurationOverlapException(box);
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
        if (t < 0) {
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
     * @return true if the Integrator samples according to a specified temperature
     */
    public boolean isIsothermal() {
        return isothermal;
    }

    /**
     * @param b specifies whether the Integrator should (if true) sample according to a specified temperature
     */
    public void setIsothermal(boolean b) {
        isothermal = b;
    }

    /**
     * @return the Box instance on which this integrator acts
     */
    public Box getBox() {
        return box;
    }

    public void saveState(Writer fw) throws IOException {
        super.saveState(fw);
        fw.write(""+currentPotentialEnergy+"\n");
    }

    public void restoreState(BufferedReader br) throws IOException {
        super.restoreState(br);
        currentPotentialEnergy = Double.parseDouble(br.readLine());
    }
}
