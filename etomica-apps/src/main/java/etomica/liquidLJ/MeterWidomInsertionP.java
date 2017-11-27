/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.action.BoxInflate;
import etomica.action.MoleculeActionTranslateTo;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorBox;
import etomica.molecule.IMolecule;
import etomica.space.Space;
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
 * 
 * The actual chemical potential can be calculated as -kT ln(<x>) where x is
 * the value returned by getDataAsScalar.
 * 
 * @author David Kofke
 */
public class MeterWidomInsertionP extends DataSourceScalar {

    public MeterWidomInsertionP(Space space, IRandom random) {
        super("exp(-\u03BC/kT)", Null.DIMENSION);//"\u03BC" is Unicode for greek "mu"
        setNInsert(100);
        atomTranslator = new MoleculeActionTranslateTo(space);
        positionSource = new RandomPositionSourceRectangular(space, random);
        boxInflate = new BoxInflate(space);
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
        energyMeter.setTarget((IMolecule)null);
        double u0 = energyMeter.getDataAsScalar();
        double v0 = box.getBoundary().volume();
        int n0 = box.getMoleculeList().getMoleculeCount();
        boxInflate.setTargetDensity(n0*n0/(v0*(n0+1)));
        boxInflate.actionPerformed();
        double u1 = energyMeter.getDataAsScalar();
        energyMeter.setTarget(testMolecule);
        if (integrator != null) temperature = integrator.getTemperature();
        for (int i = nInsert; i > 0; i--) { //perform nInsert insertions
            atomTranslator.setDestination(positionSource.randomPosition());
            atomTranslator.actionPerformed(testMolecule);
            box.addMolecule(testMolecule);
            double u = energyMeter.getDataAsScalar();
            sum += Math.exp(-epsFactor * (u+u1-u0) / temperature);
            if (Double.isInfinite(sum)) {
                throw new RuntimeException("oops");
            }
            box.removeMolecule(testMolecule);
        }
        boxInflate.setTargetDensity(n0/v0);
        boxInflate.actionPerformed();

        // ((N+1)/N)^N ~= e, so rho*e^(Z-1)
        return Math.pow((n0+1.0)/n0, n0)*v0/n0*Math.exp(-pressure*v0/n0)*sum / nInsert; //return average
    }

    /**
     * Returns the integrator associated with this class.  The box, potentialMaster
     * and temperature are taken from the integrator.
     */
    public IntegratorBox getIntegrator() {
        return integrator;
    }

    /**
     * Sets the integrator associated with this class.  The box, potentialMaster
     * and temperature are taken from the integrator.  Alternatively, you can
     * set the temperature, box and energy meter separately.
     */
    public void setIntegrator(IntegratorBox newIntegrator) {
        integrator = newIntegrator;
        setEnergyMeter(new MeterPotentialEnergy(integrator.getPotentialMaster()));
        setBox(integrator.getBox());
    }

    public void setEnergyMeter(MeterPotentialEnergy newEnergyMeter) {
        energyMeter = newEnergyMeter;
    }

    public void setBox(Box newBox) {
        this.box = newBox;
        energyMeter.setBox(box);
        positionSource.setBox(box);
        boxInflate.setBox(box);
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
        else if (integrator != null) {
            positionSource.setBox(integrator.getBox());
        }
    }
    
    /**
     * Returns the RandomPositionSource used by this meter.
     */
    public RandomPositionSource getPositionSource() {
        return positionSource;
    }

    private IntegratorBox integrator;

    /**
     * Number of insertions attempted in each call to currentValue Default is
     * 100
     */
    private int nInsert;
    private ISpecies species;
    private IMolecule testMolecule;// prototype insertion molecule
    private MoleculeActionTranslateTo atomTranslator;
    protected RandomPositionSource positionSource;
    private MeterPotentialEnergy energyMeter;
    protected Box box;
    protected double temperature;
    protected double pressure = Double.NaN;
    public double epsFactor = 1.0;
    protected BoxInflate boxInflate;
}
