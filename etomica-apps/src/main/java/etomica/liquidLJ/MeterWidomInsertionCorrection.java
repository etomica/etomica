package etomica.liquidLJ;

import etomica.action.MoleculeActionTranslateTo;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.math.function.IFunction;
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
public class MeterWidomInsertionCorrection implements IDataSource {

    public MeterWidomInsertionCorrection(Space space, IRandom random) {
        tag = new DataTag();
        dataInfo = new DataInfoDoubleArray("exp(-\u03BC/kT)", Null.DIMENSION, new int[]{4});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(4);
        setNInsert(100);
        atomTranslator = new MoleculeActionTranslateTo(space);
        positionSource = new RandomPositionSourceRectangular(space, random);
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
    
    public void setEnergyFac(double newEnergyFac) {
        uFac = newEnergyFac;
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
    
    public void setResidual(boolean newResidual) {
        residual = newResidual;
    }

    public void setPressure(double newPressureFast, double newPressureFull, IFunction vBias) {
        pressureFast = newPressureFast;
        pressureFull = newPressureFull;
        v0 = box.getBoundary().volume();
        this.vBias = vBias;
    }

    /**
     * Performs a Widom insertion average, doing nInsert insertion attempts
     * Temperature used to get exp(-uTest/kT) is that of the integrator for the
     * box
     * 
     * @return the sum of exp(-uTest/kT)/nInsert 
     */
    public IData getData() {
        energyMeterFast.setTarget((IMolecule)null);
        energyMeterFull.setTarget((IMolecule)null);
        double u0Fast = energyMeterFast.getDataAsScalar() + uFac;
        double u0Full = energyMeterFull.getDataAsScalar();
//        System.out.println(Math.exp(-(u0Full-u0Fast) / temperature));
        double sumFast = 0.0; //sum for local insertion average
        double sumFull = 0.0;
        energyMeterFast.setTarget(testMolecule);
        energyMeterFull.setTarget(testMolecule);
        double[] x = data.getData();
        for (int i = nInsert; i > 0; i--) { //perform nInsert insertions
            atomTranslator.setDestination(positionSource.randomPosition());
            atomTranslator.actionPerformed(testMolecule);
            box.addMolecule(testMolecule);
            double uFast = epsFactor*energyMeterFast.getDataAsScalar();
            sumFast += Math.exp(-uFast / temperature);
            double uFull = epsFactor*energyMeterFull.getDataAsScalar();
            sumFull += Math.exp(-uFull / temperature);;
            box.removeMolecule(testMolecule);
        }
        
        double dx = u0Full-u0Fast;
        double fac = 1;
        if (!residual) {
            // multiply by V/N
            fac = box.getBoundary().volume() / (box.getNMolecules(species)+1);
        }
        else if (!Double.isNaN(pressureFull)) {
            //  ??????
            fac = pressureFull*box.getBoundary().volume() / ((box.getNMolecules(species) + 1)*temperature);
            throw new RuntimeException("you get to figure this out.  or include ideal gas part");
        }
        double wfac = 1;
        if (!Double.isNaN(pressureFull)) {
            dx += (pressureFull-pressureFast)*(box.getBoundary().volume() - v0);
            wfac = vBias.f(box.getBoundary().volume());
        }

        x[2] = Math.exp(-dx / temperature)/wfac;
        x[0] = fac*x[2]*sumFull / nInsert;
        x[1] = fac*sumFast / nInsert;
        x[3] = fac*sumFull / nInsert;
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public void setEnergyMeter(MeterPotentialEnergy newEnergyMeterFast, MeterPotentialEnergy newEnergyMeterFull) {
        energyMeterFast = newEnergyMeterFast;
        energyMeterFull = newEnergyMeterFull;
    }
    
    public void setBox(Box newBox) {
        this.box = newBox;
        energyMeterFast.setBox(box);
        energyMeterFull.setBox(box);
        positionSource.setBox(box);
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
    protected final DataInfoDoubleArray dataInfo;
    protected final DataDoubleArray data;
    protected final DataTag tag;
    private int nInsert;
    private ISpecies species;
    private IMolecule testMolecule;// prototype insertion molecule
    private MoleculeActionTranslateTo atomTranslator;
    protected RandomPositionSource positionSource;
    private MeterPotentialEnergy energyMeterFast, energyMeterFull;
    protected Box box;
    protected double temperature;
    protected double uFac = 0;
    public double epsFactor = 1.0;
    protected boolean residual = false;
    protected double pressureFast = Double.NaN, pressureFull = Double.NaN;
    protected double v0 = -1;
    protected IFunction vBias;
}
