/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.action.MoleculeActionTranslateTo;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.species.ISpecies;
import etomica.units.dimensions.Null;
import etomica.util.random.IRandom;

/**
 * Meter to measure the chemical potential (as its exponent: exp(-mu/kT)) of a
 * species via the Widom insertion method. Call to getDataAsScalar returns a
 * Widom-insertion average (i.e., sum of exp(-energy/kT) using a configurable
 * number of insertion trials) for a molecule of the desired species in the
 * given box in its current configuration. It is expected that repeated calls
 * will be averaged to give an overall insertion average for many configurations
 * of the box.
 * <br>
 * Meter can be configured to measure the residual chemical potential -- above
 * the ideal-gas value -- or the full chemical potential.  If not configured for
 * the residual chemical potential, the insertion average will be multiplied by
 * the species number density in the box in its current configuration.  Default
 * behavior will give residual chemical potential.
 * <p>
 * The actual chemical potential can be calculated as -kT ln(<x>) where x is
 * the value returned by getDataAsScalar.
 *
 * @author David Kofke
 */
public class MeterWidomInsertionFasterer extends DataSourceScalar {

    public MeterWidomInsertionFasterer(Box box, IRandom random, PotentialCompute potentialCompute, double temperature) {
        super("exp(-\u03BC/kT)", Null.DIMENSION);//"\u03BC" is Unicode for greek "mu"
        setNInsert(100);
        setResidual(true);
        atomTranslator = new MoleculeActionTranslateTo(box.getSpace());
        positionSource = new RandomPositionSourceRectangular(box.getSpace(), random);
        this.box = box;
        positionSource.setBox(box);
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
     * Sets the species, takes a prototype molecule, and gets handle to
     * appropriate species agent in box
     */
    public void setSpecies(ISpecies s) {
        species = s;
        testMolecule = s.makeMolecule();
    }

    /**
     * Accessor for the species for which chemical potential is evaluated
     */
    public ISpecies getSpecies() {
        return species;
    }

    /**
     * Number of Widom insertions attempted with each call to currentValue
     */
    public void setNInsert(int n) {
        nInsert = n;
    }

    /**
     * Accessor to number of Widom insertions attempted with each call to
     * currentValue
     */
    public int getNInsert() {
        return nInsert;
    }

    public void setPressure(double newPressure) {
        pressure = newPressure;
    }

    /**
     * Performs a Widom insertion average, doing nInsert insertion attempts
     * Temperature used to get exp(-uTest/kT) is that of the integrator for the
     * box
     *
     * @return the sum of exp(-uTest/kT)/nInsert, multiplied by V/N if
     * <code>residual</code> is false
     */
    public double getDataAsScalar() {
        double sum = 0.0; //sum for local insertion average
        for (int i = nInsert; i > 0; i--) { //perform nInsert insertions
            atomTranslator.setDestination(positionSource.randomPosition());
            atomTranslator.actionPerformed(testMolecule);
            box.addMolecule(testMolecule);
            double u = potentialCompute.computeOneMolecule(testMolecule);
            sum += Math.exp(-u / temperature);
            if (Double.isInfinite(sum)) {
                throw new RuntimeException("oops");
            }
            box.removeMolecule(testMolecule);
        }

        if (!residual) {
            // multiply by V/N
            sum *= box.getBoundary().volume() / (box.getNMolecules(species) + 1);
        } else if (!Double.isNaN(pressure)) {
            sum *= pressure * box.getBoundary().volume() / ((box.getNMolecules(species) + 1) * temperature);
        }
        return sum / nInsert; //return average
    }

    public void setTemperature(double newTemperature) {
        this.temperature = newTemperature;
    }

    /**
     * Sets a new RandomPositionSource for this meter to use.  By default, a
     * position source is used which assumes rectangular boundaries.
     */
    public void setPositionSource(RandomPositionSource newPositionSource) {
        positionSource = newPositionSource;
        if (box != null) {
            positionSource.setBox(box);
        }
    }

    /**
     * Returns the RandomPositionSource used by this meter.
     */
    public RandomPositionSource getPositionSource() {
        return positionSource;
    }


    /**
     * Number of insertions attempted in each call to currentValue Default is
     * 100
     */
    private int nInsert;
    private ISpecies species;
    private IMolecule testMolecule;// prototype insertion molecule
    private boolean residual; // flag to specify if total or residual chemical
    // potential evaluated. Default true
    private final MoleculeActionTranslateTo atomTranslator;
    private final PotentialCompute potentialCompute;
    protected RandomPositionSource positionSource;
    protected Box box;
    protected double temperature;
    protected double pressure = Double.NaN;
}
