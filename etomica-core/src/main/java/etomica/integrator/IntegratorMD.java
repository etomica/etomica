/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.api.*;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.*;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.exception.ConfigurationOverlapException;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Time;
import etomica.util.Debug;
import etomica.util.EnumeratedType;
import etomica.util.random.IRandom;

/**
 * Superclass of all molecular-dynamics integrators.
 * Extends the Integrator class by adding methods that 
 * set the time step.
 */
public abstract class IntegratorMD extends IntegratorBox implements BoxEventListener {

    public IntegratorMD(PotentialMaster potentialMaster, IRandom random,
                        double timeStep, double temperature, Space _space) {
        super(potentialMaster,temperature);
        this.random = random;
        this.space = _space;
        setTimeStep(timeStep);
        thermostat = ThermostatType.ANDERSEN;
        setThermostatInterval(100);
        meterKE = new MeterKineticEnergy();
        atomActionRandomizeVelocity = new AtomActionRandomizeVelocity(temperature, random);
        momentum = space.makeVector();
        temperatureVec = space.makeVector();
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
    
    public void setBox(Box p) {
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

        if (thermostat == ThermostatType.HYBRID_MC) {
            oldPositionAgentManager = new AtomLeafAgentManager<Vector>(new VectorSource(space), box, Vector.class);
        }
        if (integratorMC != null) {
            integratorMC.setBox(box);
        }
    }

    protected void setup() {
        super.setup();
        currentTime = 0;
        thermostatCount = 1;
        doThermostatInternal();
    }

    public void resetStepCount() {
        super.resetStepCount();
        currentTime = 0;
    }

    /**
     * reset the integrator's kinetic energy tracker
     */
    public void reset() {
        ConfigurationOverlapException overlapException = null;
        try {
            super.reset();
        }
        catch (ConfigurationOverlapException e) {
            overlapException = e;
        }
        currentKineticEnergy = meterKE.getDataAsScalar();
        if (overlapException != null) {
            throw overlapException;
        }
        if (thermostat == ThermostatType.HYBRID_MC && isothermal && oldPositionAgentManager == null) {
            oldPositionAgentManager = new AtomLeafAgentManager<Vector>(new VectorSource(space), box, Vector.class);
        }
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
            doThermostatInternal();
        }
        if (!isothermal && integratorMC != null) {
            throw new RuntimeException("still have an integratorMC!");
        }
    }
    
    public void setTemperature(double t) {
        if (t == temperature) return;
        super.setTemperature(t);
        if (initialized) {
            atomActionRandomizeVelocity.setTemperature(temperature);
            // trigger immediate thermostat
            thermostatCount = 1;
            doThermostatInternal();
        }

        if (integratorMC != null) {
            integratorMC.setTemperature(temperature);
        }
    }
    
    public static class ThermostatType extends EnumeratedType {
        protected ThermostatType(String label) {super(label);}       
        private static final long serialVersionUID = 1L;
        public static final ThermostatType VELOCITY_SCALING = new ThermostatType("Velocity Scaling");
        public static final ThermostatType ANDERSEN = new ThermostatType("Anderson");
        public static final ThermostatType ANDERSEN_SINGLE = new ThermostatType("Andersen Single");
        public static final ThermostatType ANDERSEN_SCALING = new ThermostatType("Andersen Scaling");
        public static final ThermostatType HYBRID_MC = new ThermostatType("Hybrid MC");
        //public static final ThermostatType NOSE_HOOVER;
        public static ThermostatType[] choices() {
            return new ThermostatType[] {VELOCITY_SCALING,ANDERSEN,ANDERSEN_SINGLE,HYBRID_MC,ANDERSEN_SCALING};
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

    public void setIntegratorMC(IntegratorMC integratorMC, int mcSteps) {
        if (thermostat != ThermostatType.HYBRID_MC || !isothermal) {
            throw new RuntimeException("integratorMC only works with HYBRID MC thermostat");
        }
        this.integratorMC = integratorMC;
        integratorMC.setTemperature(temperature);
        if (box != null){
            integratorMC.setBox(box);
            integratorMC.reset();       	
        }
        this.mcSteps = mcSteps;
    }

    public IntegratorMC getIntegratorMC() {
        return integratorMC;
    }

    /**
     * Sets the type of thermostat used by the integrator.
     * @param aThermostat the desired thermostat
     */
    public void setThermostat(ThermostatType aThermostat) {
        thermostat = aThermostat;
        if (thermostat == ThermostatType.HYBRID_MC && box != null) {
            oldPositionAgentManager = new AtomLeafAgentManager<Vector>(new VectorSource(space), box, Vector.class);
        }
        else if (thermostat != ThermostatType.HYBRID_MC) {
            if (oldPositionAgentManager != null) {
                oldPositionAgentManager.dispose();
                oldPositionAgentManager = null;
            }
            if (integratorMC != null) {
                throw new RuntimeException("still have an integratorMC!");
            }
        }
        if (thermostat == ThermostatType.ANDERSEN_SCALING) {
            setThermostatNoDrift(true);
        }
        oldEnergy = Double.NaN;
    }

    /**
     * Configures the thermostat to remove any net momenta imparted
     * to the system.
     */
    public void setThermostatNoDrift(boolean newNoDrift) {
        thermostatNoDrift = newNoDrift;
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

    public int getThermostatInterval() {
        return thermostatInterval;
    }

    /**
     * Fires the thermostat if the appropriate interval has been reached.
     */
    protected void doThermostatInternal() {
        if (--thermostatCount == 0) {
            doThermostat();
        }
    }

    /**
     * thermostat implementation.  This method takes whatever action is appropriate
     * for the integrator's thermostat and updates the state of the integrator.
     */
    public void doThermostat() {
        thermostatting = true;
        thermostatCount = thermostatInterval;
        boolean rejected = false;
        if (thermostat == ThermostatType.ANDERSEN || thermostat == ThermostatType.ANDERSEN_SCALING || thermostat == ThermostatType.HYBRID_MC || !initialized) {
            if (thermostat == ThermostatType.HYBRID_MC) {
                if (!Double.isNaN(oldEnergy)) {
                    // decide whether or not to go back to the old configuration
                    double newPotentialEnergy = meterPE.getDataAsScalar();
                    double newKineticEnergy = meterKE.getDataAsScalar();
//                    System.out.println(newPotentialEnergy+" "+newKineticEnergy+" "+oldEnergy);
                    double energyDiff = newPotentialEnergy + newKineticEnergy - oldEnergy;
//                    System.out.println(energyDiff+" "+Math.exp(-energyDiff/temperature));
                    if (energyDiff > 0 && Math.exp(-energyDiff/temperature) < random.nextDouble()) {
                        // energy increased and we are rejecting the trajectory
                        rejected = true;
                        IAtomList leafAtoms = box.getLeafList();
                        for (int i=0; i<leafAtoms.getAtomCount(); i++) {
                            IAtom a = leafAtoms.getAtom(i);
                            a.getPosition().E(oldPositionAgentManager.getAgent(a));
                        }
                        oldEnergy = oldPotentialEnergy;
//                        System.out.println("rejected "+energyDiff+" => "+oldEnergy);
//                        System.out.println(" *** check *** "+meterPE.getDataAsScalar());
                        nRejected++;
                    }
                    else {
                        // accepting the trajectory.  save positions
                        IAtomList leafAtoms = box.getLeafList();
                        for (int i=0; i<leafAtoms.getAtomCount(); i++) {
                            IAtom a = leafAtoms.getAtom(i);
                            oldPositionAgentManager.getAgent(a).E(a.getPosition());
                        }
                        oldPotentialEnergy = newPotentialEnergy;
                        oldEnergy = newPotentialEnergy;
//                        System.out.println("accepted "+energyDiff+" => "+oldEnergy);

                        nAccepted++;
                    }
                }
                else {
                    // initialize with the current configuration
                    IAtomList leafAtoms = box.getLeafList();
                    for (int i=0; i<leafAtoms.getAtomCount(); i++) {
                        IAtom a = leafAtoms.getAtom(i);
                        oldPositionAgentManager.getAgent(a).E(a.getPosition());
                    }
                    oldPotentialEnergy = meterPE.getDataAsScalar();
                    oldEnergy = oldPotentialEnergy;
                }

                if (initialized && integratorMC != null) {
                    for (int i=0; i<mcSteps; i++) {
                        integratorMC.doStep();
                    }
                    IAtomList leafAtoms = box.getLeafList();
                    for (int i=0; i<leafAtoms.getAtomCount(); i++) {
                        IAtom a = leafAtoms.getAtom(i);
                        oldPositionAgentManager.getAgent(a).E(a.getPosition());
                    }
                    randomizeMomenta();
                    if (thermostatNoDrift) {
                        shiftMomenta();
                    }
                    reset();
                    oldPotentialEnergy = currentPotentialEnergy;
                    oldEnergy = oldPotentialEnergy;
                }
                else if (rejected) {
                    // we've put the atoms back, need to reset (force recalc)
                    randomizeMomenta();
                    if (thermostatNoDrift) {
                        shiftMomenta();
                    }
                    reset();
                }
                else {
                    randomizeMomenta();
                    if (thermostatNoDrift) {
                        shiftMomenta();
                    }
                    currentKineticEnergy = meterKE.getDataAsScalar();
                }
//                System.out.print(oldEnergy+" ");
                oldEnergy += currentKineticEnergy;
//                System.out.println(" ===> "+oldEnergy);
            }
            else if (thermostat == ThermostatType.ANDERSEN || !initialized) {
                randomizeMomenta();
                if (thermostatNoDrift) {
                    shiftMomenta();
                }
                currentKineticEnergy = meterKE.getDataAsScalar();
            }
            else if (thermostat == ThermostatType.ANDERSEN_SCALING) {
                randomizeTotalKE();
                currentKineticEnergy = meterKE.getDataAsScalar();
            }
        }
        if (thermostat == ThermostatType.VELOCITY_SCALING || !isothermal) {
            shiftMomenta();
            scaleMomenta();
        }
        else if (thermostat == ThermostatType.ANDERSEN_SINGLE) {
            if (initialized) {
                IAtomList atomList = box.getLeafList();
                int atomCount = atomList.getAtomCount();
                if (atomCount > 0) {
                    int index = random.nextInt(atomList.getAtomCount());
                    IAtomKinetic a = (IAtomKinetic)atomList.getAtom(index);
                    double m = a.getType().getMass();
                    if (m == Double.POSITIVE_INFINITY) return;
                    currentKineticEnergy -= 0.5*m*a.getVelocity().squared();
                    randomizeMomentum(a);
                    currentKineticEnergy += 0.5*m*a.getVelocity().squared();
                }
            }
        }
        // ANDERSEN was handled at the start
        else if (thermostat != ThermostatType.ANDERSEN && thermostat != ThermostatType.HYBRID_MC && thermostat != ThermostatType.ANDERSEN_SCALING) {
            throw new RuntimeException("Unknown thermostat: "+thermostat);
        }
        thermostatting = false;
    }

    public double getHybridAcceptance() {
        return ((double)nAccepted)/(nAccepted+nRejected);
    }

    protected void randomizeTotalKE() {
        shiftMomenta();

        IAtomList leafList= box.getLeafList();
        int nLeaf = box.getLeafList().getAtomCount();
        int nSkipped = 0;
        double totalMass = 0;
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic atom = (IAtomKinetic)leafList.getAtom(iLeaf);
            double mass = atom.getType().getMass();
            if(mass == Double.POSITIVE_INFINITY || mass == 0) nSkipped--;
            else totalMass += mass;
        }
        nLeaf -= nSkipped;
        for (int i=0; i<space.D(); i++) {
            double sumBKE = 0;
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                double x = random.nextGaussian();
                sumBKE += x*x;
            }
            temperatureVec.setX(i, sumBKE);
        }
        temperatureVec.TE(temperature/nLeaf);
        scaleMomenta(temperatureVec);
        
        if (!thermostatNoDrift && totalMass > 0) {
            // pick net velocity from Maxwell-Boltzmann distribution
            momentum.E(0);
            for (int i=0; i<space.D(); i++) {
                double x = random.nextGaussian();
                double BKE = x*x;
                double v = 2*Math.sqrt(BKE*temperature/totalMass);
                momentum.setX(i, v);
            }
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtomKinetic atom = (IAtomKinetic)leafList.getAtom(iLeaf);
                atom.getVelocity().PE(momentum);
            }
        }
    }
    
    /**
     * randomizes the velocities for the given box using velocities
     * chosen form a Maxwell-Boltzmann distribution as in the Andersen 
     * thermostat.  The state of the integrator needs to be updated 
     * after calling this method.
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
        atomActionRandomizeVelocity.actionPerformed(atom);
    }
    
    /**
     * Subtracts velocity from all atoms such that the new net momentum of
     * the whole system is 0.  This could also be implemented to subtract
     * momentum from each atom, but then 0-mass atoms would explode.
     */
    protected void shiftMomenta() {
        momentum.E(0);
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        if (nLeaf == 0) return;
        if (nLeaf > 1) {
            double totalMass = 0;
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtom a = leafList.getAtom(iLeaf);
                double mass = a.getType().getMass();
                if (mass != Double.POSITIVE_INFINITY) {
                    momentum.PEa1Tv1(mass,((IAtomKinetic)a).getVelocity());
                    totalMass += mass;
                }
            }
            if (totalMass == 0) return;
            momentum.TE(1.0/totalMass);
            //momentum is now net velocity
            //set net momentum to 0
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
                double rm = a.getType().rm();
                if (rm != 0 && rm != Double.POSITIVE_INFINITY) {
                    a.getVelocity().ME(momentum);
                }
            }
            if (Debug.ON) {
                momentum.E(0);
                for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
                    double mass = a.getType().getMass();
                    if (mass != Double.POSITIVE_INFINITY) {
                        momentum.PEa1Tv1(mass,a.getVelocity());
                    }
                }
                momentum.TE(1.0/totalMass);
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
        temperatureVec.E(temperature);
        scaleMomenta(temperatureVec);
    }

    protected void scaleMomenta(Vector t) {
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        currentKineticEnergy = 0;
        if (nLeaf == 0) return;
        // calculate current kinetic temperature.
        for (int i = 0; i < space.D(); i++) {
            // scale independently in each dimension
            double sum = 0.0;
            int nLeafNotFixed = 0;
            for (int iAtom = 0; iAtom<nLeaf; iAtom++) {
                IAtomKinetic atom = (IAtomKinetic)leafList.getAtom(iAtom);
                double mass = atom.getType().getMass();
                if(mass == Double.POSITIVE_INFINITY) continue;
                double v = atom.getVelocity().getX(i);
                sum += mass*v*v;
                nLeafNotFixed++;
            }
            if (sum == 0) {
                if (t.getX(i) == 0) {
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
            double s = Math.sqrt(t.getX(i) / (sum / nLeafNotFixed));
            currentKineticEnergy += 0.5*sum*s*s;
            if (s == 1) continue;
            for (int iAtom = 0; iAtom<nLeaf; iAtom++) {
                IAtomKinetic atom = (IAtomKinetic)leafList.getAtom(iAtom);
                Vector vel = atom.getVelocity();
                vel.setX(i, vel.getX(i)*s); //scale momentum
            }
        }
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
    
    public void boxMoleculeAdded(BoxMoleculeEvent e) {

        IMolecule mole = e.getMolecule();
        IAtomList atomList = mole.getChildList();
        for(int i = 0; i < atomList.getAtomCount(); i++) {
            IAtom a = atomList.getAtom(i);
            if (a instanceof IAtomKinetic) {
                randomizeMomentum((IAtomKinetic)a);
            }
        }

    }
    
    public void boxMoleculeRemoved(BoxMoleculeEvent e) { }

    public void boxGlobalAtomLeafIndexChanged(BoxIndexEvent e) { }
    
    public void boxAtomLeafIndexChanged(BoxAtomIndexEvent e) { }
    
    public void boxMoleculeIndexChanged(BoxMoleculeIndexEvent e) { }
    
    public void boxNumberMolecules(BoxMoleculeCountEvent e) { }
    
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
    protected final Vector momentum;
    protected double currentTime;
    protected final Space space;
    protected boolean thermostatting = false;
    protected boolean thermostatNoDrift = false;

    protected double oldEnergy = Double.NaN, oldPotentialEnergy = Double.NaN;
    protected long nRejected = 0, nAccepted = 0;
    protected AtomLeafAgentManager<Vector> oldPositionAgentManager = null;

    protected IntegratorMC integratorMC;
    protected int mcSteps;
    protected final Vector temperatureVec;

    public static class VectorSource implements AgentSource<Vector> {

        protected final Space space;
        
        public VectorSource(Space space) {
            this.space = space;
        }

        public Vector makeAgent(IAtom a, Box agentBox) {
            Vector p = space.makeVector();
            p.E(a.getPosition());
            return p;
        }

        public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {}
    }
}

