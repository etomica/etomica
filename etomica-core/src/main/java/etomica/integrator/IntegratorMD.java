/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.action.AtomActionRandomizeVelocity;
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
import etomica.meta.annotations.IgnoreProperty;
import etomica.molecule.IMolecule;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Dimensioned;
import etomica.units.dimensions.Time;
import etomica.util.Debug;
import etomica.util.EnumeratedType;
import etomica.util.random.IRandom;

/**
 * Superclass of all molecular-dynamics integrators.
 */
public abstract class IntegratorMD extends IntegratorBox implements BoxEventListener {

    protected final IRandom random;
    protected final Vector momentum;
    protected final Space space;
    protected final Vector temperatureVec;
    protected double timeStep;
    protected double currentKineticEnergy;
    protected ThermostatType thermostat;
    protected int thermostatCount, thermostatInterval;
    protected DataSourceScalar meterKE;
    protected AtomActionRandomizeVelocity atomActionRandomizeVelocity;
    protected MeterTemperature meterTemperature;
    protected double currentTime;
    protected boolean thermostatting = false;
    protected boolean thermostatNoDrift = false;
    protected double oldEnergy = Double.NaN, oldPotentialEnergy = Double.NaN;
    protected long nRejected = 0, nAccepted = 0;
    protected AtomLeafAgentManager<Vector> oldPositionAgentManager = null;
    protected IntegratorMC integratorMC;
    protected int mcSteps;

    /**
     * Constructs integrator with a default for non-isothermal sampling.
     *
     * @param potentialMaster PotentialMaster instance used to compute the energy and forces
     * @param random          random number generator used for initial velocities and some thermostats
     * @param timeStep        time step for integration
     * @param temperature     used by thermostat and/or to initialize velocities
     * @param space           the governing space, used to generate vectors
     */
    public IntegratorMD(PotentialMaster potentialMaster, IRandom random,
                        double timeStep, double temperature, Space space) {
        super(potentialMaster, temperature);
        this.random = random;
        this.space = space;
        setTimeStep(timeStep);
        thermostat = ThermostatType.ANDERSEN;
        setThermostatInterval(100);
        meterKE = new MeterKineticEnergy();
        atomActionRandomizeVelocity = new AtomActionRandomizeVelocity(temperature, random);
        momentum = this.space.makeVector();
        temperatureVec = this.space.makeVector();
    }

    /**
     * @return the integration time step
     */
    @Dimensioned(dimension = Time.class)
    public final double getTimeStep() {
        return timeStep;
    }

    /**
     * @param t the new integration time step
     */
    public void setTimeStep(double t) {
        timeStep = t;
    }


    public void setBox(Box box) {
        if (this.box != null) {
            this.box.getEventManager().removeListener(this);
        }
        super.setBox(box);
        if (meterKE instanceof MeterKineticEnergy) {
            ((MeterKineticEnergy) meterKE).setBox(box);
        }
        meterTemperature = new MeterTemperature(box, space.D());
        meterTemperature.setKineticEnergyMeter(meterKE);
        this.box.getEventManager().addListener(this);

        if (thermostat == ThermostatType.HYBRID_MC) {
            oldPositionAgentManager = new AtomLeafAgentManager<Vector>(new VectorSource(space), this.box, Vector.class);
        }
        if (integratorMC != null) {
            integratorMC.setBox(this.box);
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

    public void reset() {
        ConfigurationOverlapException overlapException = null;
        try {
            super.reset();
        } catch (ConfigurationOverlapException e) {
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

    protected void doStepInternal() {
        currentTime += timeStep;
    }

    /**
     * @return the amount of simulation time that has elapsed since construction or since last
     * call to resetStepCount
     */
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

    /**
     * A Monte Carlo integrator is employed with this MD integrator when performing hybrid MD/MC
     * sampling. The HYBRIC_MC thermostat will accept/reject a segment of the MD trajectory, and
     * then, optionally, the IntegratorMC will perform additional sampling (e.g., grand-canonical),
     * leaving the system in a state ready for the next MD segment.
     *
     * @param integratorMC the Monte Carlo integrator to be used for hybrid MC/MD
     * @param mcSteps      number of Monte Carlo trial steps to perform after MD segment
     * @throws RuntimeException if HYBRID_MC thermostat is not enabled or if isothermal is false
     */
    public void setIntegratorMC(IntegratorMC integratorMC, int mcSteps) {
        if (thermostat != ThermostatType.HYBRID_MC || !isothermal) {
            throw new RuntimeException("integratorMC only works with HYBRID MC thermostat");
        }
        this.integratorMC = integratorMC;
        integratorMC.setTemperature(temperature);
        if (box != null) {
            integratorMC.setBox(box);
            integratorMC.reset();
        }
        this.mcSteps = mcSteps;
    }

    /**
     * @return the Monte Carlo integrator used to enable hybrid MD/MC sampling. Will be null
     * if not previously set via setIntegratorMC.
     */
    public IntegratorMC getIntegratorMC() {
        return integratorMC;
    }

    /**
     * Configures the thermostat to remove any net momenta imparted
     * to the system after use of the thermostat.
     *
     * @param newNoDrift if true, will remove net momentum after thermostat use
     */
    public void setThermostatNoDrift(boolean newNoDrift) {
        thermostatNoDrift = newNoDrift;
    }

    /**
     *
     * @return flag specifying whether net momentum should be removed after application of thermostat
     */
    public boolean isThermostatNoDrift() {return thermostatNoDrift;}

    /**
     *
     * @return the type of thermostat used to implement isothermal sampling.
     */
    public ThermostatType getThermostat() {
        return thermostat;
    }

    /**
     * Sets the type of thermostat used by the integrator for isothermal sampling. The default
     * at construction is ANDERSEN.
     *
     * @param aThermostat the desired thermostat
     */
    public void setThermostat(ThermostatType aThermostat) {
        thermostat = aThermostat;
        if (thermostat == ThermostatType.HYBRID_MC && box != null) {
            oldPositionAgentManager = new AtomLeafAgentManager<Vector>(new VectorSource(space), box, Vector.class);
        } else if (thermostat != ThermostatType.HYBRID_MC) {
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
     *
     * @return the value of the thermostat interval
     */
    public int getThermostatInterval() {
        return thermostatInterval;
    }

    /**
     * Sets the number of integrator steps between thermostat
     * actions (for thermostats such as velocity scaling and Andersen).
     *
     * @param interval number of integrator steps between thermostat
     *                 activity
     */
    public void setThermostatInterval(int interval) {
        if (interval < 1) throw new IllegalArgumentException("Thermostat Interval must be positive");
        thermostatInterval = interval;
        thermostatCount = interval;
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
     * Invokes thermostat, modifying velocities according to specified temperature and
     * algorithm of thermostat specified via setThermostat.
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
                    if (energyDiff > 0 && Math.exp(-energyDiff / temperature) < random.nextDouble()) {
                        // energy increased and we are rejecting the trajectory
                        rejected = true;
                        IAtomList leafAtoms = box.getLeafList();
                        for (int i = 0; i < leafAtoms.getAtomCount(); i++) {
                            IAtom a = leafAtoms.getAtom(i);
                            a.getPosition().E(oldPositionAgentManager.getAgent(a));
                        }
                        oldEnergy = oldPotentialEnergy;
//                        System.out.println("rejected "+energyDiff+" => "+oldEnergy);
//                        System.out.println(" *** check *** "+meterPE.getDataAsScalar());
                        nRejected++;
                    } else {
                        // accepting the trajectory.  save positions
                        IAtomList leafAtoms = box.getLeafList();
                        for (int i = 0; i < leafAtoms.getAtomCount(); i++) {
                            IAtom a = leafAtoms.getAtom(i);
                            oldPositionAgentManager.getAgent(a).E(a.getPosition());
                        }
                        oldPotentialEnergy = newPotentialEnergy;
                        oldEnergy = newPotentialEnergy;
//                        System.out.println("accepted "+energyDiff+" => "+oldEnergy);

                        nAccepted++;
                    }
                } else {
                    // initialize with the current configuration
                    IAtomList leafAtoms = box.getLeafList();
                    for (int i = 0; i < leafAtoms.getAtomCount(); i++) {
                        IAtom a = leafAtoms.getAtom(i);
                        oldPositionAgentManager.getAgent(a).E(a.getPosition());
                    }
                    oldPotentialEnergy = meterPE.getDataAsScalar();
                    oldEnergy = oldPotentialEnergy;
                }

                if (initialized && integratorMC != null) {
                    for (int i = 0; i < mcSteps; i++) {
                        integratorMC.doStep();
                    }
                    IAtomList leafAtoms = box.getLeafList();
                    for (int i = 0; i < leafAtoms.getAtomCount(); i++) {
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
                } else if (rejected) {
                    // we've put the atoms back, need to reset (force recalc)
                    randomizeMomenta();
                    if (thermostatNoDrift) {
                        shiftMomenta();
                    }
                    reset();
                } else {
                    randomizeMomenta();
                    if (thermostatNoDrift) {
                        shiftMomenta();
                    }
                    currentKineticEnergy = meterKE.getDataAsScalar();
                }
//                System.out.print(oldEnergy+" ");
                oldEnergy += currentKineticEnergy;
//                System.out.println(" ===> "+oldEnergy);
            } else if (thermostat == ThermostatType.ANDERSEN || !initialized) {
                randomizeMomenta();
                if (thermostatNoDrift) {
                    shiftMomenta();
                }
                currentKineticEnergy = meterKE.getDataAsScalar();
            } else if (thermostat == ThermostatType.ANDERSEN_SCALING) {
                randomizeTotalKE();
                currentKineticEnergy = meterKE.getDataAsScalar();
            }
        }
        if (thermostat == ThermostatType.VELOCITY_SCALING || !isothermal) {
            shiftMomenta();
            scaleMomenta();
        } else if (thermostat == ThermostatType.ANDERSEN_SINGLE) {
            if (initialized) {
                IAtomList atomList = box.getLeafList();
                int atomCount = atomList.getAtomCount();
                if (atomCount > 0) {
                    int index = random.nextInt(atomList.getAtomCount());
                    IAtomKinetic a = (IAtomKinetic) atomList.getAtom(index);
                    double m = a.getType().getMass();
                    if (m == Double.POSITIVE_INFINITY) return;
                    currentKineticEnergy -= 0.5 * m * a.getVelocity().squared();
                    randomizeMomentum(a);
                    currentKineticEnergy += 0.5 * m * a.getVelocity().squared();
                }
            }
        }
        // ANDERSEN was handled at the start
        else if (thermostat != ThermostatType.ANDERSEN && thermostat != ThermostatType.HYBRID_MC && thermostat != ThermostatType.ANDERSEN_SCALING) {
            throw new RuntimeException("Unknown thermostat: " + thermostat);
        }
        thermostatting = false;
    }

    /**
     *
     * @return fraction of Monte Carlo trials that were accepted when using HYBRID_MC thermostat,
     * since construction, or since last call to resetHybridAcceptance.
     */
    public double getHybridAcceptance() {
        return ((double) nAccepted) / (nAccepted + nRejected);
    }

    /**
     * Zeros counters used to track acceptance fraction of HYBRID_MC thermostat.
     */
    public void resetHybridAcceptance() {
        nAccepted = 0;
        nRejected = 0;
    }

    /**
     * Scales all velocities so that the total kinetic energy matches a value selected at random for the
     * previously specified temperature.
     */
    protected void randomizeTotalKE() {
        shiftMomenta();

        IAtomList leafList = box.getLeafList();
        int nLeaf = box.getLeafList().getAtomCount();
        int nSkipped = 0;
        double totalMass = 0;
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic atom = (IAtomKinetic) leafList.getAtom(iLeaf);
            double mass = atom.getType().getMass();
            if (mass == Double.POSITIVE_INFINITY || mass == 0) nSkipped--;
            else totalMass += mass;
        }
        nLeaf -= nSkipped;
        for (int i = 0; i < space.D(); i++) {
            double sumBKE = 0;
            for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
                double x = random.nextGaussian();
                sumBKE += x * x;
            }
            temperatureVec.setX(i, sumBKE);
        }
        temperatureVec.TE(temperature / nLeaf);
        scaleMomenta(temperatureVec);

        if (!thermostatNoDrift && totalMass > 0) {
            // pick net velocity from Maxwell-Boltzmann distribution
            momentum.E(0);
            for (int i = 0; i < space.D(); i++) {
                double x = random.nextGaussian();
                double BKE = x * x;
                double v = 2 * Math.sqrt(BKE * temperature / totalMass);
                momentum.setX(i, v);
            }
            for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
                IAtomKinetic atom = (IAtomKinetic) leafList.getAtom(iLeaf);
                atom.getVelocity().PE(momentum);
            }
        }
    }

    /**
     * Randomizes the velocities of all atoms in the box using values
     * chosen from a Maxwell-Boltzmann distribution, as in the Andersen
     * thermostat. Does not reset or otherwise updated integrator for
     * new velocities (some integrators, e.g. IntegratorHard, will need to be
     * reset after calling this method).
     */
    protected void randomizeMomenta() {
        atomActionRandomizeVelocity.setTemperature(temperature);
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            atomActionRandomizeVelocity.actionPerformed(leafList.getAtom(iLeaf));
        }
    }

    /**
     * Randomizes the velocities of the given atom using value
     * chosen from a Maxwell-Boltzmann distribution, as in the Andersen
     * thermostat. Does not reset or otherwise updated integrator for
     * new velocities (some integrators, e.g. IntegratorHard, will need to be
     * reset after calling this method).
     *
     * @param atom whose momenta is be randomized
     *
     */
    protected void randomizeMomentum(IAtomKinetic atom) {
        atomActionRandomizeVelocity.setTemperature(temperature);
        atomActionRandomizeVelocity.actionPerformed(atom);
    }

    /**
     * Subtracts velocity from all atoms such that the new net momentum of
     * the whole system is zero.
     */
    protected void shiftMomenta() {
        momentum.E(0);
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        if (nLeaf == 0) return;
        if (nLeaf > 1) {
            double totalMass = 0;
            for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
                IAtom a = leafList.getAtom(iLeaf);
                double mass = a.getType().getMass();
                if (mass != Double.POSITIVE_INFINITY) {
                    momentum.PEa1Tv1(mass, ((IAtomKinetic) a).getVelocity());
                    totalMass += mass;
                }
            }
            if (totalMass == 0) return;
            momentum.TE(1.0 / totalMass);
            //momentum is now net velocity
            //set net momentum to 0
            for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
                IAtomKinetic a = (IAtomKinetic) leafList.getAtom(iLeaf);
                double rm = a.getType().rm();
                if (rm != 0 && rm != Double.POSITIVE_INFINITY) {
                    a.getVelocity().ME(momentum);
                }
            }
            if (Debug.ON) {
                momentum.E(0);
                for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic) leafList.getAtom(iLeaf);
                    double mass = a.getType().getMass();
                    if (mass != Double.POSITIVE_INFINITY) {
                        momentum.PEa1Tv1(mass, a.getVelocity());
                    }
                }
                momentum.TE(1.0 / totalMass);
                if (Math.sqrt(momentum.squared()) > 1.e-10) {
                    System.out.println("Net momentum per leaf atom is " + momentum + " but I expected it to be 0");
                }
            }
            momentum.E(0);
        }
    }

    /**
     * Crude method to enforce constant-temperature constraint, by
     * scaling momenta of all atoms by a constant factor so that
     * box adheres to setpoint temperature.  The state of the
     * integrator may need to be updated after calling this method.
     */
    protected void scaleMomenta() {
        temperatureVec.E(temperature);
        scaleMomenta(temperatureVec);
    }

    // used internally
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
            for (int iAtom = 0; iAtom < nLeaf; iAtom++) {
                IAtomKinetic atom = (IAtomKinetic) leafList.getAtom(iAtom);
                double mass = atom.getType().getMass();
                if (mass == Double.POSITIVE_INFINITY) continue;
                double v = atom.getVelocity().getX(i);
                sum += mass * v * v;
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
                    throw new RuntimeException("atoms have no velocity component in " + i + " dimension");
                }
                randomizeMomenta();
                i--;
                // try again, we could infinite loop in theory
                continue;
            }
            double s = Math.sqrt(t.getX(i) / (sum / nLeafNotFixed));
            currentKineticEnergy += 0.5 * sum * s * s;
            if (s == 1) continue;
            for (int iAtom = 0; iAtom < nLeaf; iAtom++) {
                IAtomKinetic atom = (IAtomKinetic) leafList.getAtom(iAtom);
                Vector vel = atom.getVelocity();
                vel.setX(i, vel.getX(i) * s); //scale momentum
            }
        }
    }

    /**
     * @return the temperature meter used for velocity rescaling.
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
        for (int i = 0; i < atomList.getAtomCount(); i++) {
            IAtom a = atomList.getAtom(i);
            if (a instanceof IAtomKinetic) {
                randomizeMomentum((IAtomKinetic) a);
            }
        }

    }

    public void boxMoleculeRemoved(BoxMoleculeEvent e) {
    }

    public void boxGlobalAtomLeafIndexChanged(BoxIndexEvent e) {
    }

    public void boxAtomLeafIndexChanged(BoxAtomIndexEvent e) {
    }

    public void boxMoleculeIndexChanged(BoxMoleculeIndexEvent e) {
    }

    public void boxNumberMolecules(BoxMoleculeCountEvent e) {
    }

    public static class ThermostatType extends EnumeratedType {
        public static final ThermostatType VELOCITY_SCALING = new ThermostatType("Velocity Scaling");
        public static final ThermostatType ANDERSEN = new ThermostatType("Anderson");
        public static final ThermostatType ANDERSEN_SINGLE = new ThermostatType("Andersen Single");
        public static final ThermostatType ANDERSEN_SCALING = new ThermostatType("Andersen Scaling");
        public static final ThermostatType HYBRID_MC = new ThermostatType("Hybrid MC");

        protected ThermostatType(String label) {
            super(label);
        }

        //public static final ThermostatType NOSE_HOOVER;
        public static ThermostatType[] choices() {
            return new ThermostatType[]{VELOCITY_SCALING, ANDERSEN, ANDERSEN_SINGLE, HYBRID_MC, ANDERSEN_SCALING};
        }
    }

    /**
     * An AgentSource that permits another vector to be associated with each atom.
     * This is needed by most MD integrators save information needed to advance
     * the atom positions/velocities.
     */
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

        public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {
        }
    }
}

