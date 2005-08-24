package etomica.integrator;

import etomica.Simulation;
import etomica.action.AtomActionRandomizeVelocity;
import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.exception.MethodNotImplementedException;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.ICoordinateKinetic;
import etomica.units.Dimension;
import etomica.util.Constants;
import etomica.util.Default;
import etomica.util.Constants.TypedConstant;
/**
 * Superclass of all molecular-dynamics integrators.
 * Extends the Integrator class by adding methods that 
 * set the time step.
 */

public abstract class IntegratorMD extends Integrator {
    
    public IntegratorMD(PotentialMaster potentialMaster) {
        super(potentialMaster);
        setTimeStep(Default.TIME_STEP);
        thermostat = ANDERSEN;
        atomIterator = new AtomIteratorLeafAtoms();
        setThermostatInterval(100);
        meterKE = new MeterKineticEnergy();
        atomActionRandomizeVelocity = new AtomActionRandomizeVelocity();
        atomActionRandomizeVelocity.setTemperature(temperature);
        meterTemperature = new MeterTemperature();
        currentKineticEnergy = new double[1];
    }

    /**
     * Sets integration time step.
     * Updates zero-point counters used internally to manage the elapsedTime method
     */
    public void setTimeStep(double t) {
        timeStep = t;
    }
    public final double getTimeStep() {return timeStep;}
    public Dimension getTimeStepDimension() {return Dimension.TIME;}

    protected void setup() {
        super.setup();
        thermostatCount = 1;
        meterKE.setPhase(firstPhase);
        doThermostat();
    }
    
    /**
     * reset the integrator's kinetic energy tracker
     */
    public void reset() {
        super.reset();
        meterKE.setPhase(firstPhase);
        if (phase.length != currentKineticEnergy.length) {
            currentKineticEnergy = new double[phase.length];
        }
        for (int i=0; i<phase.length; i++) {
            meterKE.setPhase(phase[i]);
            currentKineticEnergy[i] = meterKE.getDataAsScalar();
        }
    }

    /**
     * @return the current kinetic energy as tracked by the integrator
     */
    public double[] getKineticEnergy() {
        return currentKineticEnergy;
    }

    public void setIsothermal(boolean b) {
        super.setIsothermal(b);
        if (initialized && isothermal) {
            // trigger immediate thermostat
            thermostatCount = 1;
            doThermostat();
        }
    }
    
    public void setTemperature(double t) {
        if (t == temperature) return;
        super.setTemperature(t);
        if (initialized) {
            atomActionRandomizeVelocity.setTemperature(temperature);
            // trigger immediate thermostat
            thermostatCount = 1;
            doThermostat();
        }
    }
    
    public static class ThermostatType extends TypedConstant {
        protected ThermostatType(String label) {super(label);}       
        public Constants.TypedConstant[] choices() {return CHOICES;}
    }
    protected static final ThermostatType[] CHOICES = 
        new ThermostatType[] {
            new ThermostatType("Velocity Scaling"), new ThermostatType("Andersen"),
            new ThermostatType("Andersen Single"), new ThermostatType("Nose Hoover")};
    public static final ThermostatType VELOCITY_SCALING = CHOICES[0];
    public static final ThermostatType ANDERSEN = CHOICES[1];
    public static final ThermostatType ANDERSEN_SINGLE = CHOICES[2];
    public static final ThermostatType NOSE_HOOVER = CHOICES[3];
    
    /**
     * Sets the type of thermostat used by the integrator.
     * @param aThermostat the desired thermostat
     */
    public void setThermostat(ThermostatType aThermostat) {
        thermostat = aThermostat;
    }

    /**
     * Sets the number of integrator intervals between thermostat
     * actions (for velocity scaling and Andersen thermostat).
     * @param interval number of integrator intervals between thermostat
     * activity
     */
    public void setThermostatInterval(int interval) {
        if (interval < 1) throw new IllegalArgumentException("Thermostat Interval must be positive");
        thermostatInterval = interval;
        thermostatCount = interval;
    }

    /**
     * thermostat implementation.  This method takes whatever action is appropriate
     * for the integrator's thermostat and updates the state of the integrator.
     */
    public void doThermostat() {
        if (--thermostatCount == 0) {
            thermostatCount = thermostatInterval;
            if (thermostat == VELOCITY_SCALING || !isothermal) {
                // calculate current kinetic temperature
                for (int i=0; i<phase.length; i++) {
                    scaleMomenta(phase[i]);
                    meterKE.setPhase(phase[i]);
                    currentKineticEnergy[i] = meterKE.getDataAsScalar();
                }
            }
            else if (thermostat == ANDERSEN) {
                for (int i=0; i<phase.length; i++) {
                    randomizeMomenta(phase[i]);
                    meterKE.setPhase(phase[i]);
                    currentKineticEnergy[i] = meterKE.getDataAsScalar();
                }
            }
            else if (thermostat == ANDERSEN_SINGLE) {
                for (int i=0; i<phase.length; i++) {
                    AtomList atomList = phase[i].getSpeciesMaster().atomList;
                    int index = Simulation.random.nextInt(atomList.size());
                    Atom a = atomList.get(index);
                    double m = ((AtomTypeLeaf)a.type).getMass();
                    currentKineticEnergy[i] -= 0.5*m*((ICoordinateKinetic)a.coord).velocity().squared();
                    randomizeMomentum(atomList.get(index));
                    currentKineticEnergy[i] += 0.5*m*((ICoordinateKinetic)a.coord).velocity().squared();
                }
            }
            else if (thermostat == NOSE_HOOVER) {
                throw new MethodNotImplementedException("feel free to write the Nose Hoover thermostat");
            }
            else {
                throw new RuntimeException("Unknown thermostat: "+thermostat);
            }
        }
    }
    
    /**
     * randomizes the velocities for the given phase using velocities
     * chosen form a Maxwell-Boltzmann distribution as in the Andersen 
     * thermostat.  The state of the integrator needs to be updated 
     * after calling this method.
     * @param aPhase phase whose atomic momenta are to be randomized
     */
    protected void randomizeMomenta(Phase aPhase) {
        atomIterator.setPhase(aPhase);
        atomActionRandomizeVelocity.setTemperature(temperature);
        atomIterator.allAtoms(atomActionRandomizeVelocity);
    }
    
    /**
     * randomizes the velocity of an atom in the given phase using velocities
     * chosen form a Maxwell-Boltzmann distribution as in the Andersen 
     * thermostat.  The state of the integrator needs to be updated 
     * after calling this method.
     * @param atom whose momenta is be randomized
     */
    protected void randomizeMomentum(Atom atom) {
        atomActionRandomizeVelocity.setTemperature(temperature);
        atomActionRandomizeVelocity.actionPerformed(atom);
    }
    
    /**
     * Crude method to enforce constant-temperature constraint
     * Scales momenta of all atoms by a constant factor so that 
     * phase adheres to setpoint temperature.  The state of the 
     * integrator needs to be updated after calling this method.
     * @return the factor velocities were scaled by 
     */
    protected double scaleMomenta(Phase aPhase) {
        atomIterator.setPhase(aPhase);
        atomIterator.reset();
        // calculate current kinetic temperature
        meterTemperature.setPhase(aPhase);
        double t = meterTemperature.getDataAsScalar();
        if (t == temperature) return 1.0;
        double s = Math.sqrt(temperature / t);
        double scale = s;
        if (t == 0) {
            randomizeMomenta(aPhase);
            meterTemperature.setPhase(aPhase);
            t = meterTemperature.getDataAsScalar();
            s = Math.sqrt(temperature / t);
        }
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            ((ICoordinateKinetic)a.coord).velocity().TE(s); //scale momentum
        }
        return scale;
    }

    /**
     * returns the temperature meter used for velocity rescaling.
     */
    public DataSourceScalar getMeterTemperature() {
        return meterTemperature;
    }
    /**
     * Sets the temperature meter used to calculate temperature for
     * velocity rescaling.  You only need to call this method if 
     * the standard MeterTemperature won't work.
     */
    public void setMeterTemperature(MeterTemperature meter) {
        meterTemperature = meter;
    }
    
    /**
     * Elementary time step for the MD simulation
     */
    protected double timeStep;
    protected double[] currentKineticEnergy;
    protected AtomIteratorLeafAtoms atomIterator;
    protected ThermostatType thermostat;
    private int thermostatCount, thermostatInterval;
    protected MeterKineticEnergy meterKE;
    private AtomActionRandomizeVelocity atomActionRandomizeVelocity;
    private MeterTemperature meterTemperature;
    
}//end of IntegratorMD
    