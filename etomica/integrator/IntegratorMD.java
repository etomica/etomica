package etomica.integrator;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.exception.ConfigurationOverlapException;
import etomica.exception.MethodNotImplementedException;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateKinetic;
import etomica.units.Dimension;
import etomica.util.EnumeratedType;
/**
 * Superclass of all molecular-dynamics integrators.
 * Extends the Integrator class by adding methods that 
 * set the time step.
 */

public abstract class IntegratorMD extends IntegratorPhase {
    
    public IntegratorMD(PotentialMaster potentialMaster, double timeStep, double temperature) {
        super(potentialMaster,temperature);
        setTimeStep(timeStep);
        thermostat = ANDERSEN;
        atomIterator = new AtomIteratorLeafAtoms();
        setThermostatInterval(100);
        meterKE = new MeterKineticEnergy();
        atomActionRandomizeVelocity = new AtomActionRandomizeVelocity(temperature);
        meterTemperature = new MeterTemperature();
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
        try {
            super.setup();
        }
        catch (ConfigurationOverlapException e) {}
        thermostatCount = 1;
        meterKE.setPhase(phase);
        doThermostat();
    }
    
    /**
     * reset the integrator's kinetic energy tracker
     */
    public void reset() throws ConfigurationOverlapException {
        meterKE.setPhase(phase);
        currentKineticEnergy = meterKE.getDataAsScalar();
        super.reset();
    }

    /**
     * @return the current kinetic energy as tracked by the integrator
     */
    public double getKineticEnergy() {
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
    
    public static class ThermostatType extends EnumeratedType {
        protected ThermostatType(String label) {super(label);}       
        public EnumeratedType[] choices() {return CHOICES;}

        /**
         * Required to guarantee singleton when deserializing.
         * @return the singleton INSTANCE
         */
        private Object readResolve() {
            for (int i=0; i<CHOICES.length; i++) {
                if (this.toString().equals(CHOICES[i].toString())) {
                    return CHOICES[i];
                }
            }
            throw new RuntimeException("unknown thermostat type: "+this);
        }
    }
    protected static final ThermostatType[] CHOICES = 
        new ThermostatType[] {
            new ThermostatType("Velocity Scaling"), new ThermostatType("Andersen"),
            new ThermostatType("Andersen Single"), new ThermostatType("Nose Hoover (unimplemented)")};
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

    public ThermostatType getThermostat() {
        return thermostat;
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
            if (thermostat == ANDERSEN || !initialized) {
                // if initializing the system always randomize the velocity
                randomizeMomenta(phase);
                meterKE.setPhase(phase);
                currentKineticEnergy = meterKE.getDataAsScalar();
            }
            if (thermostat == VELOCITY_SCALING || (!initialized && thermostat == NOSE_HOOVER)) {
                // rescale randomized velocities for Nose Hoover during initialization
                scaleMomenta(phase);
                meterKE.setPhase(phase);
                currentKineticEnergy = meterKE.getDataAsScalar();
            }
            else if (thermostat == ANDERSEN_SINGLE) {
                AtomList atomList = phase.getSpeciesMaster().atomList;
                int index = Simulation.random.nextInt(atomList.size());
                Atom a = atomList.get(index);
                double m = ((AtomTypeLeaf)a.type).getMass();
                currentKineticEnergy -= 0.5*m*((ICoordinateKinetic)a.coord).velocity().squared();
                randomizeMomentum(atomList.get(index));
                currentKineticEnergy += 0.5*m*((ICoordinateKinetic)a.coord).velocity().squared();
            }
            else if (thermostat == NOSE_HOOVER) {
                throw new MethodNotImplementedException("feel free to write the Nose Hoover thermostat");
            }
            // ANDERSEN was handled at the start
            else if (thermostat != ANDERSEN) {
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
    protected double currentKineticEnergy;
    protected AtomIteratorLeafAtoms atomIterator;
    protected ThermostatType thermostat;
    private int thermostatCount, thermostatInterval;
    protected MeterKineticEnergy meterKE;
    private AtomActionRandomizeVelocity atomActionRandomizeVelocity;
    private MeterTemperature meterTemperature;
    
}//end of IntegratorMD
    