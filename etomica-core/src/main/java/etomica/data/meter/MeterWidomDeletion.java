/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.compute.PotentialCompute;
import etomica.species.ISpecies;
import etomica.units.dimensions.Null;
import etomica.util.random.IRandom;

/**
 * Meter to measure the chemical potential (as its exponent: exp(-mu/kT)) of a
 * species via the Widom deletion method. Call to getDataAsScalar returns a
 * Widom-deletion average (i.e., sum of exp(-energy/kT) using a configurable
 * number of deletion trials) for a molecule of the desired species in the
 * given box in its current configuration. It is expected that repeated calls
 * will be averaged to give an overall deletion average for many configurations
 * of the box.
 * <br>
 * Meter can be configured to measure the residual chemical potential -- above
 * the ideal-gas value -- or the full chemical potential. If not configured for
 * the residual chemical potential, the deletion average will be multiplied by
 * the species number density in the box in its current configuration. Default
 * behavior will give residual chemical potential.
 * <p>
 * The actual chemical potential can be calculated as -kT ln(<x>) where x is
 * the value returned by getDataAsScalar.
 *
 * @author Based on MeterWidomInsertion by David Kofke
 */
public class MeterWidomDeletion extends DataSourceScalar {

    public MeterWidomDeletion(Box box, IRandom random, PotentialCompute potentialCompute, double temperature) {
        super("exp(-\u03BC/kT)", Null.DIMENSION);//"\u03BC" is Unicode for greek "mu"
        setNDelete(100);
        setResidual(true);
        this.box = box;
        this.random = random;
        this.temperature = temperature;
        this.potentialCompute = potentialCompute;
    }

    /**
     * Sets flag specifying if full or residual chemical potential is computed
     * Default is <code>true</code> (only residual is computed)
     */
    public void setResidual(boolean b) {
        residual = b;
    }

    /**
     * Accessor for flag specifying if full or residual chemical potential is
     * computed
     */
    public boolean isResidual() {
        return residual;
    }

    /**
     * Sets the species and gets handle to appropriate species agent in box
     */
    public void setSpecies(ISpecies s) {
        species = s;
    }

    /**
     * Accessor for the species for which chemical potential is evaluated
     */
    public ISpecies getSpecies() {
        return species;
    }

    /**
     * Number of Widom deletions attempted with each call to currentValue
     */
    public void setNDelete(int n) {
        nDelete = n;
    }

    /**
     * Accessor to number of Widom deletions attempted with each call to
     * currentValue
     */
    public int getNDelete() {
        return nDelete;
    }

    public void setPressure(double newPressure) {
        pressure = newPressure;
    }

    /**
     * Performs a Widom deletion average, doing nDelete deletion attempts
     * Temperature used to get exp(-uTest/kT) is that specified for this meter
     *
     * @return the sum of exp(-uTest/kT)/nDelete, multiplied by N/V if
     * <code>residual</code> is false
     */
    public double getDataAsScalar() {
        double sum = 0.0; //sum for local deletion average
        int nMolecules = box.getNMolecules(species);

        if (nMolecules == 0) {
            return Double.NaN; // Cannot perform deletions if there are no molecules
        }

        for (int i = nDelete; i > 0; i--) { //perform nDelete deletion tests
            // Randomly select a molecule from the box
            int randomIndex = random.nextInt(nMolecules);
            IMoleculeList moleculeList = box.getMoleculeList(species);
            IMolecule selectedMolecule = moleculeList.get(randomIndex);

            // Calculate the energy of this molecule in the current configuration
            double u = potentialCompute.computeOneMolecule(selectedMolecule);

            sum += Math.exp(+u / temperature);
            if (Double.isInfinite(sum)) {
                throw new RuntimeException("oops");
            }
        }

        if (!residual) {
            // multiply by N/V for full chemical potential (note: using N instead of N+1)
            sum *= box.getNMolecules(species) / box.getBoundary().volume();
        } else if (!Double.isNaN(pressure)) {
            // For NPT simulations with specified pressure
            sum *= box.getNMolecules(species) * temperature / (pressure * box.getBoundary().volume());
        }

        return sum / nDelete; //return average
    }

    public void setTemperature(double newTemperature) {
        this.temperature = newTemperature;
    }

    /**
     * Number of deletions attempted in each call to currentValue Default is
     * 100
     */
    private int nDelete;
    private ISpecies species;
    private boolean residual; // flag to specify if total or residual chemical
    // potential evaluated. Default true
    private final PotentialCompute potentialCompute;
    protected Box box;
    protected IRandom random;
    protected double temperature;
    protected double pressure = Double.NaN;
}