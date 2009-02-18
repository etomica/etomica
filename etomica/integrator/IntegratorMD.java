package etomica.integrator;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.api.IAtomKinetic;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IBoxMoleculeAddedEvent;
import etomica.api.IEvent;
import etomica.api.IListener;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.box.BoxMoleculeAddedEvent;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.space.ISpace;
import etomica.units.Dimension;
import etomica.units.Time;
import etomica.util.Debug;
import etomica.util.EnumeratedType;

/**
 * Superclass of all molecular-dynamics integrators.
 * Extends the Integrator class by adding methods that 
 * set the time step.
 */
public abstract class IntegratorMD extends IntegratorBox implements IListener {

    public IntegratorMD(IPotentialMaster potentialMaster, IRandom random, 
            double timeStep, double temperature, ISpace _space) {
        super(potentialMaster,temperature);
        this.random = random;
        this.space = _space;
        setTimeStep(timeStep);
        thermostat = ThermostatType.ANDERSEN;
        setThermostatInterval(100);
        meterKE = new MeterKineticEnergy();
        atomActionRandomizeVelocity = new AtomActionRandomizeVelocity(temperature, random);
        momentum = space.makeVector();
    }

    /**
     * Sets integration time step.
     * Updates zero-point counters used internally to manage the elapsedTime method
     */
    public void setTimeStep(double t) {
        timeStep = t;
    }
    public final double getTimeStep() {return timeStep;}
    public Dimension getTimeStepDimension() {return Time.DIMENSION;}
    
    public void setBox(IBox p) {
        if (box != null) {
            box.getEventManager().removeListener(this);
        }
        super.setBox(p);
        if (meterKE instanceof MeterKineticEnergy) {
            ((MeterKineticEnergy)meterKE).setBox(p);
        }
        meterTemperature = new MeterTemperature(p, space.D());
        meterTemperature.setKineticEnergyMeter(meterKE);
        box.getEventManager().addListener(this);
    }

    protected void setup() {
        super.setup();
        currentTime = 0;
        thermostatCount = 1;
        doThermostat();
    }

    public void resetStepCount() {
        super.resetStepCount();
        currentTime = 0;
    }

    /**
     * reset the integrator's kinetic energy tracker
     */
    public void reset() {
        super.reset();
        currentKineticEnergy = meterKE.getDataAsScalar();
    }

    public void doStepInternal() {
        currentTime += timeStep;
    }
    
    public double getCurrentTime() {
        return currentTime;
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
        private static final long serialVersionUID = 1L;
        public static final ThermostatType VELOCITY_SCALING = new ThermostatType("Velocity Scaling");
        public static final ThermostatType ANDERSEN = new ThermostatType("Anderson");
        public static final ThermostatType ANDERSEN_NODRIFT = new ThermostatType("Anderson without drift");
        public static final ThermostatType ANDERSEN_SINGLE = new ThermostatType("Andersen Single");
        //public static final ThermostatType NOSE_HOOVER;
        public static ThermostatType[] choices() {
            return new ThermostatType[] {VELOCITY_SCALING,ANDERSEN,ANDERSEN_SINGLE,ANDERSEN_NODRIFT};
        }

        /**
         * Required to guarantee singleton when deserializing.
         * @return the singleton INSTANCE
         */
        private Object readResolve() {
            ThermostatType[] choices = choices();
            for (int i=0; i<choices.length; i++) {
                if (this.toString().equals(choices[i].toString())) {
                    return choices[i];
                }
            }
            throw new RuntimeException("unknown thermostat type: "+this);
        }
    }
    
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
            if (thermostat == ThermostatType.ANDERSEN || thermostat == ThermostatType.ANDERSEN_NODRIFT || !initialized) {
                // if initializing the system always randomize the velocity
                randomizeMomenta();
                if (thermostat == ThermostatType.ANDERSEN) {
                    currentKineticEnergy = meterKE.getDataAsScalar();
                }
            }
            if (thermostat == ThermostatType.VELOCITY_SCALING || thermostat == ThermostatType.ANDERSEN_NODRIFT || !isothermal) {
                shiftMomenta();
                if (thermostat == ThermostatType.VELOCITY_SCALING || !isothermal) {
                    scaleMomenta();
                }
                else {
                    // shifting changes the kinetic energy, so we need to recalculate it
                    // (scaleMomenta does this on its own)
                    currentKineticEnergy = meterKE.getDataAsScalar();
                }
            }
            else if (thermostat == ThermostatType.ANDERSEN_SINGLE) {
                if (initialized) {
                    IAtomList atomList = box.getLeafList();
                    int atomCount = atomList.getAtomCount();
                    if (atomCount > 0) {
                        int index = random.nextInt(atomList.getAtomCount());
                        IAtomKinetic a = (IAtomKinetic)atomList.getAtom(index);
                        double m = ((IAtom)a).getType().getMass();
                        currentKineticEnergy -= 0.5*m*a.getVelocity().squared();
                        randomizeMomentum(a);
                        currentKineticEnergy += 0.5*m*a.getVelocity().squared();
                    }
                }
            }
            // ANDERSEN was handled at the start
            else if (thermostat != ThermostatType.ANDERSEN) {
                throw new RuntimeException("Unknown thermostat: "+thermostat);
            }
        }
    }
    
    /**
     * randomizes the velocities for the given box using velocities
     * chosen form a Maxwell-Boltzmann distribution as in the Andersen 
     * thermostat.  The state of the integrator needs to be updated 
     * after calling this method.
     * @param aBox box whose atomic momenta are to be randomized
     */
    protected void randomizeMomenta() {
        atomActionRandomizeVelocity.setTemperature(temperature);
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            atomActionRandomizeVelocity.actionPerformed(leafList.getAtom(iLeaf));
        }
    }
    
    /**
     * randomizes the velocity of an atom in the given box using velocities
     * chosen form a Maxwell-Boltzmann distribution as in the Andersen 
     * thermostat.  The state of the integrator needs to be updated 
     * after calling this method.
     * @param atom whose momenta is be randomized
     */
    protected void randomizeMomentum(IAtomKinetic atom) {
        atomActionRandomizeVelocity.setTemperature(temperature);
        atomActionRandomizeVelocity.actionPerformed((IAtom)atom);
    }
    
    /**
     * Subtracts net momentum from all atoms such that the new net momentum of
     * the whole system is 0.
     */
    protected void shiftMomenta() {
        momentum.E(0);
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        if (nLeaf == 0) return;
        if (nLeaf > 1) {
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtom a = leafList.getAtom(iLeaf);
                double mass = a.getType().getMass();
                if (mass != Double.POSITIVE_INFINITY) {
                    momentum.PEa1Tv1(mass,((IAtomKinetic)a).getVelocity());
                }
            }
            momentum.TE(1.0/nLeaf);
            //set net momentum to 0
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
                double rm = ((IAtom)a).getType().rm();
                if (rm != 0) {
                    a.getVelocity().PEa1Tv1(-rm,momentum);
                }
            }
            if (Debug.ON) {
                momentum.E(0);
                for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
                    double mass = ((IAtom)a).getType().getMass();
                    if (mass != Double.POSITIVE_INFINITY) {
                        momentum.PEa1Tv1(mass,a.getVelocity());
                    }
                }
                momentum.TE(1.0/nLeaf);
                if (Math.sqrt(momentum.squared()) > 1.e-10) {
                    System.out.println("Net momentum per leaf atom is "+momentum+" but I expected it to be 0");
                }
            }
            momentum.E(0);
        }
    }

    /**
     * Crude method to enforce constant-temperature constraint
     * Scales momenta of all atoms by a constant factor so that 
     * box adheres to setpoint temperature.  The state of the 
     * integrator may need to be updated after calling this method.
     */
    protected void scaleMomenta() {
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        currentKineticEnergy = 0;
        if (nLeaf == 0) return;
        // calculate current kinetic temperature.
        for (int i = 0; i < space.D(); i++) {
            // scale independently in each dimension
            double sum = 0.0;
            for (int iAtom = 0; iAtom<nLeaf; iAtom++) {
                IAtomKinetic atom = (IAtomKinetic)leafList.getAtom(iAtom);
                double mass = ((IAtom)atom).getType().getMass();
                if(mass == Double.POSITIVE_INFINITY) continue;
                double v = atom.getVelocity().x(i);
                sum += mass*v*v;
            }
            if (sum == 0) {
                if (temperature == 0) {
                    continue;
                }
                if (i > 0) {
                    // wonky.  possible in theory.  but then, you called
                    // scaleMomenta, so you're probably a bad person and
                    // deserve this.
                    throw new RuntimeException("atoms have no velocity component in "+i+" dimension");
                }
                randomizeMomenta();
                i--;
                // try again, we could infinite loop in theory
                continue;
            }
            double s = Math.sqrt(temperature / (sum / nLeaf));
            currentKineticEnergy += 0.5*sum*s*s;
            if (s == 1) continue;
            for (int iAtom = 0; iAtom<nLeaf; iAtom++) {
                IAtomKinetic atom = (IAtomKinetic)leafList.getAtom(iAtom);
                IVectorMutable vel = atom.getVelocity(); 
                vel.setX(i, vel.x(i)*s); //scale momentum
            }
        }
        
        
//        double t = meterTemperature.getDataAsScalar();
//        if (t == temperature) return;
//        double s = Math.sqrt(temperature / t);
//        double scale = s;
//        if (t == 0) {
//            randomizeMomenta();
//            t = meterTemperature.getDataAsScalar();
//            s = Math.sqrt(temperature / t);
//        }
//        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
//            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
//            a.getVelocity().TE(s); //scale momentum
//        }
//        return;
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
    
    public void actionPerformed(IEvent event) {
        if (event instanceof IBoxMoleculeAddedEvent) {
            IMolecule mole = ((IBoxMoleculeAddedEvent)event).getMolecule();
            IAtomList atomList = mole.getChildList();
            for(int i = 0; i < atomList.getAtomCount(); i++) {
                randomizeMomentum((IAtomKinetic)atomList.getAtom(i));
            }
        }
    }

    private static final long serialVersionUID = 2L;
    /**
     * Elementary time step for the MD simulation
     */
    protected final IRandom random;
    protected double timeStep;
    protected double currentKineticEnergy;
    protected ThermostatType thermostat;
    protected int thermostatCount, thermostatInterval;
    protected DataSourceScalar meterKE;
    protected AtomActionRandomizeVelocity atomActionRandomizeVelocity;
    protected MeterTemperature meterTemperature;
    protected final IVectorMutable momentum;
    protected double currentTime;
    protected final ISpace space;
}

