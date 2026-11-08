/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.heatcapacity;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.integrator.ActionZeroMomentum;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.IMolecule;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2Harmonic;
import etomica.potential.P2WCA;
import etomica.potential.P3BondAngle;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

import java.util.Arrays;
import java.util.Collections;

/**
 * Dilute gas of monatomic, diatomic and bent triatomic molecules, for measuring the heat capacity
 * and comparing it with the equipartition theorem.  Atoms interact via WCA between molecules, and
 * molecules are flexible (harmonic bonds and a harmonic bend).  Reduced units:
 * sigma = epsilon = atom mass = kB = 1.
 */
public class HeatCapacity extends Simulation {

    public static final double BOX_LENGTH = 15;
    /** shorter than the atom diameter, so that bonded atoms visibly overlap */
    public static final double BOND_LENGTH = 0.7;
    public static final double BOND_ANGLE = Math.PI * 105.0 / 180.0;
    /** gives the same bend frequency as stiffness 50 with bond length 1 */
    public static final double ANGLE_STIFFNESS = 25;
    public static final int DEFAULT_MOLECULES = 100;
    /** largest time step, set by the WCA collisions rather than by any intramolecular mode */
    public static final double MAX_TIME_STEP = 0.005;

    public final SpeciesGeneral speciesMono, speciesDi, speciesTri;
    public final Box box;
    public final PotentialMasterBonding pcBonding;
    public final PotentialMasterList pcPair;
    public final PotentialComputeAggregate pcAggregate;
    public final IntegratorVelocityVerlet integrator;
    public final ActivityIntegrate activityIntegrate;
    public final ConfigurationLattice config;
    public final HeatCapacityMeter meter;

    protected final P2Harmonic p2Bond;
    protected final P3BondAngle p3Bend;
    protected double bondStiffness;

    public HeatCapacity() {
        this(1.0, 200, 0, DEFAULT_MOLECULES, 0);
    }

    public HeatCapacity(double temperature, double bondStiffness, int nMono, int nDi, int nTri) {
        super(Space3D.getInstance());

        AtomType typeMono = AtomType.simple("M", 1.0);
        AtomType typeDi = AtomType.simple("D", 1.0);
        AtomType typeTri = AtomType.simple("T", 1.0);

        speciesMono = SpeciesGeneral.monatomic(space, typeMono, true);
        double halfAngle = 0.5 * BOND_ANGLE;
        speciesDi = new SpeciesBuilder(space)
                .addAtom(typeDi, Vector.of(-0.5 * BOND_LENGTH, 0, 0))
                .addAtom(typeDi, Vector.of(0.5 * BOND_LENGTH, 0, 0))
                .setDynamic(true)
                .build();
        speciesTri = new SpeciesBuilder(space)
                .addAtom(typeTri, Vector.of(-BOND_LENGTH * Math.sin(halfAngle), BOND_LENGTH * Math.cos(halfAngle), 0))
                .addAtom(typeTri, Vector.of(0, 0, 0))
                .addAtom(typeTri, Vector.of(BOND_LENGTH * Math.sin(halfAngle), BOND_LENGTH * Math.cos(halfAngle), 0))
                .setDynamic(true)
                .build();
        addSpecies(speciesMono);
        addSpecies(speciesDi);
        addSpecies(speciesTri);

        box = makeBox(new BoundaryRectangularPeriodic(space, BOX_LENGTH));
        box.setNMolecules(speciesMono, nMono);
        box.setNMolecules(speciesDi, nDi);
        box.setNMolecules(speciesTri, nTri);

        pcBonding = new PotentialMasterBonding(getSpeciesManager(), box);
        this.bondStiffness = bondStiffness;
        p2Bond = new P2Harmonic(bondStiffness, BOND_LENGTH);
        p3Bend = new P3BondAngle(BOND_ANGLE, ANGLE_STIFFNESS);
        pcBonding.setBondingPotentialPair(speciesDi, p2Bond, Collections.singletonList(new int[]{0, 1}));
        pcBonding.setBondingPotentialPair(speciesTri, p2Bond, Arrays.asList(new int[]{0, 1}, new int[]{1, 2}));
        pcBonding.setBondingPotentialTriplet(speciesTri, p3Bend, Collections.singletonList(new int[]{0, 1, 2}));

        // the bonding info excludes all intramolecular pairs, including the 1-3 pair of the triatomic
        pcPair = new PotentialMasterList(getSpeciesManager(), box, 2, 2.5, pcBonding.getBondingInfo());
        P2WCA p2WCA = new P2WCA(1.0, 1.0);
        AtomType[] types = {typeMono, typeDi, typeTri};
        for (int i = 0; i < types.length; i++) {
            for (int j = i; j < types.length; j++) {
                pcPair.setPairPotential(types[i], types[j], p2WCA);
            }
        }
        pcAggregate = new PotentialComputeAggregate(pcBonding, pcPair);

        integrator = new IntegratorThermostatted(pcAggregate, timeStepFor(bondStiffness), temperature);
        integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        integrator.setThermostatNoDrift(true);
        integrator.setIsothermal(true);

        // restore each molecule's conformation before placing it, because a flexible molecule might
        // be distorted, or have its atoms wrapped to opposite sides of the box
        config = new ConfigurationLattice(new LatticeCubicFcc(space), space) {
            public void initializeCoordinates(Box b) {
                for (IMolecule molecule : b.getMoleculeList()) {
                    molecule.getType().initializeConformation(molecule);
                }
                super.initializeCoordinates(b);
            }
        };
        config.initializeCoordinates(box);

        meter = new HeatCapacityMeter(this);
        integrator.getEventManager().addListener(meter);

        activityIntegrate = new ActivityIntegrate(integrator, Long.MAX_VALUE, true);
        getController().addActivity(activityIntegrate);
    }

    /**
     * Returns the time step, 1/40 of the period of the fastest intramolecular mode, which is the
     * larger of the antisymmetric stretch and the bend of the triatomic (for equal masses,
     * w^2 = 2.26 k and w^2 = 2 * 2.26 * k_angle / l^2, respectively).  The diatomic stretch
     * (w^2 = 2 k) is always slower.  Composition is not considered, so that a given stiffness
     * gives the same step for all compositions.  The step is capped for the sake of WCA collisions.
     */
    public static double timeStepFor(double k) {
        double f = 1 - Math.cos(BOND_ANGLE) + 1;
        double w2stretch = f * k;
        double w2bend = 2 * f * ANGLE_STIFFNESS / (BOND_LENGTH * BOND_LENGTH);
        double period = 2 * Math.PI / Math.sqrt(Math.max(w2stretch, w2bend));
        return Math.min(period / 40, MAX_TIME_STEP);
    }

    public double getBondStiffness() {
        return bondStiffness;
    }

    /** Sets the bond stiffness and the corresponding time step. */
    public void setBondStiffness(double k) {
        bondStiffness = k;
        p2Bond.setSpringConstant(k);
        integrator.setTimeStep(timeStepFor(k));
        if (integrator.isInitialized()) integrator.reset();
    }

    public int getN(ISpecies species) {
        return box.getNMolecules(species);
    }

    public void setN(ISpecies species, int n) {
        box.setNMolecules(species, n);
    }

    /**
     * Adds heat Q to the system by scaling the center-of-mass velocities of the molecules, which
     * leaves rotation and vibration untouched and the net momentum at zero.
     *
     * @return false if there is no translational kinetic energy to scale (so no heat was added)
     */
    public boolean addHeatToTranslation(double Q) {
        ActionZeroMomentum.zeroMomenta(box);
        Vector vCOM = space.makeVector();
        double ke = 0;
        for (IMolecule molecule : box.getMoleculeList()) {
            centerOfMassVelocity(molecule, vCOM);
            ke += 0.5 * moleculeMass(molecule) * vCOM.squared();
        }
        if (ke <= 0) return false;
        double dScale = Math.sqrt(1 + Q / ke) - 1;
        for (IMolecule molecule : box.getMoleculeList()) {
            centerOfMassVelocity(molecule, vCOM);
            IAtomList atoms = molecule.getChildList();
            for (int i = 0; i < atoms.size(); i++) {
                ((IAtomKinetic) atoms.get(i)).getVelocity().PEa1Tv1(dScale, vCOM);
            }
        }
        return true;
    }

    private static double moleculeMass(IMolecule molecule) {
        double m = 0;
        IAtomList atoms = molecule.getChildList();
        for (int i = 0; i < atoms.size(); i++) {
            m += atoms.get(i).getType().getMass();
        }
        return m;
    }

    private static void centerOfMassVelocity(IMolecule molecule, Vector vCOM) {
        vCOM.E(0);
        IAtomList atoms = molecule.getChildList();
        for (int i = 0; i < atoms.size(); i++) {
            vCOM.PEa1Tv1(atoms.get(i).getType().getMass(), ((IAtomKinetic) atoms.get(i)).getVelocity());
        }
        vCOM.TE(1.0 / moleculeMass(molecule));
    }

    /**
     * Velocity-Verlet integrator with an Andersen-type thermostat that conserves momentum.  Each
     * firing picks a random pair of atoms and redraws their relative velocity from the Maxwell
     * distribution, leaving the velocity of the pair's center of mass unchanged.  Net momentum
     * therefore stays zero without further action, and the velocities sample the canonical
     * distribution on the zero-momentum manifold.  Instead randomizing one atom and then
     * removing the resulting net momentum leaves the temperature low by a fraction 1/N.
     */
    private class IntegratorThermostatted extends IntegratorVelocityVerlet {
        private final Vector vPair = space.makeVector();
        private final Vector vRel = space.makeVector();

        IntegratorThermostatted(PotentialCompute pc, double timeStep, double temperature) {
            super(pc, HeatCapacity.this.random, timeStep, temperature, HeatCapacity.this.box);
        }

        public void reset() {
            // molecules added since the last reset have random velocities
            if (initialized) ActionZeroMomentum.zeroMomenta(box);
            super.reset();
        }

        public void doThermostat() {
            if (!initialized || !isothermal) {
                super.doThermostat();
                return;
            }
            thermostatCount = thermostatInterval;
            IAtomList atoms = box.getLeafList();
            int n = atoms.size();
            if (n < 2) return;
            int i = random.nextInt(n);
            int j = random.nextInt(n - 1);
            if (j >= i) j++;
            IAtomKinetic a = (IAtomKinetic) atoms.get(i);
            IAtomKinetic b = (IAtomKinetic) atoms.get(j);
            double ma = a.getType().getMass(), mb = b.getType().getMass();
            vPair.Ea1Tv1(ma, a.getVelocity());
            vPair.PEa1Tv1(mb, b.getVelocity());
            vPair.TE(1.0 / (ma + mb));
            double s = Math.sqrt(temperature * (ma + mb) / (ma * mb));
            for (int d = 0; d < vRel.getD(); d++) {
                vRel.setX(d, s * random.nextGaussian());
            }
            a.getVelocity().E(vPair);
            a.getVelocity().PEa1Tv1(mb / (ma + mb), vRel);
            b.getVelocity().E(vPair);
            b.getVelocity().PEa1Tv1(-ma / (ma + mb), vRel);
        }
    }
}
