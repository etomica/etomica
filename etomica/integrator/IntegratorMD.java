package etomica.integrator;

import etomica.Atom;
import etomica.Constants;
import etomica.Default;
import etomica.Integrator;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.action.AtomActionRandomizeVelocity;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterKineticEnergy;
import etomica.exception.MethodNotImplementedException;
import etomica.space.ICoordinateKinetic;
import etomica.units.Dimension;
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

    public void initialize() {
        super.initialize();
        thermostatCount = 1;
        currentKineticEnergy = new double[phase.length];
        doThermostat();
    }
    
    /**
     * reset the integrator's kinetic energy tracker
     */
    public void reset() {
        super.reset();
        meterKE.setPhase(phase);
        currentKineticEnergy = meterKE.getData();
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
        atomActionRandomizeVelocity.setTemperature(temperature);
        if (initialized) {
            // trigger immediate thermostat
            thermostatCount = 1;
            doThermostat();
        }
    }
    
    public static class ThermostatType extends etomica.DataType {
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
                for (int i=0; i<phase.length; i++) {
                    scaleMomenta(phase[i]);
                    currentKineticEnergy[i] = temperature * (phase[i].atomCount() * phase[i].boundary().dimensions().D()) / 2.0;
                }
            }
            else if (thermostat == ANDERSEN) {
                for (int i=0; i<phase.length; i++) {
                    randomizeMomenta(phase[i]);
                }
                currentKineticEnergy = meterKE.getData();
            }
            else if (thermostat == ANDERSEN_SINGLE) {
                for (int i=0; i<phase.length; i++) {
                    AtomList atomList = phase[i].speciesMaster.atomList;
                    int index = Simulation.random.nextInt(atomList.size());
                    randomizeMomentum(atomList.get(index));
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
        double t = 2.0*meterKE.getDataAsScalar(aPhase) / (aPhase.atomCount() * aPhase.space().D());
        if (t == temperature) return 1.0;
        if (t == 0) {
            randomizeMomenta(aPhase);
            t = 2.0*meterKE.getDataAsScalar(aPhase) / (aPhase.atomCount() * aPhase.space().D());
        }
        double s = Math.sqrt(temperature / t);
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            ((ICoordinateKinetic)a.coord).velocity().TE(s); //scale momentum
        }
        return s;
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
    
}//end of IntegratorMD
    