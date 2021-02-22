/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.controller.Activity;
import etomica.action.controller.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHardFasterer;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.IntegratorMDFasterer;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.cell.NeighborCellManagerFasterer;
import etomica.nbr.list.NeighborListManagerFasterer;
import etomica.nbr.list.NeighborListManagerFastererHard;
import etomica.potential.*;
import etomica.potential.compute.NeighborManagerHard;
import etomica.potential.compute.PotentialComputePair;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

public class SimGlassFasterer extends Simulation {

    public final SpeciesGeneral speciesA, speciesB;
    public final Box box;
    public final IntegratorMDFasterer integrator;
    public final TimeSourceFasterer timeSource;

    public final MCMoveSwapFasterer swapMove;
    public final IntegratorMCFasterer integratorMC;
    public final PotentialChoice potentialChoice;
    protected P2HardGeneric p2AA;

    public enum PotentialChoice {LJ, WCA, SS, HS, WS}

    public double sigmaB;

    protected int chs;

    public SimGlassFasterer(int D, int nA, int nB, double density, double temperature, boolean doSwap, PotentialChoice pc, double tStep) {
        this(D, nA, nB, density, temperature, doSwap, pc, tStep, null);
    }

    public SimGlassFasterer(int D, int nA, int nB, double density, double temperature, boolean doSwap, PotentialChoice pc, double tStep, int[] randSeeds) {
        super(Space.getInstance(D));
        if (randSeeds != null) {
            setRandom(new RandomMersenneTwister(randSeeds));
        }
        this.potentialChoice = pc;
        //species
        speciesA = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesA);
        speciesB = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesB);

        box = this.makeBox();

        NeighborListManagerFasterer neighborManager = pc == PotentialChoice.HS
                ? new NeighborListManagerFastererHard(getSpeciesManager(), box, 2, 2.99, BondingInfo.noBonding())
                : new NeighborListManagerFasterer(getSpeciesManager(), box, 2, 2.99, BondingInfo.noBonding());
        NeighborCellManagerFasterer neighborManagerMC = new NeighborCellManagerFasterer(getSpeciesManager(), box, 2, BondingInfo.noBonding());
        PotentialComputePair potentialMaster = new PotentialComputePair(this, box, neighborManager);
        PotentialComputePair potentialMasterMC = new PotentialComputePair(this, box, neighborManagerMC);

        if (potentialChoice == PotentialChoice.LJ) { //3D KA-80-20; 2D KA-65-35
            sigmaB = 0.88;
            P2LennardJones potentialAA = new P2LennardJones(space);
            P2SoftSphericalTruncated p2TruncatedAA = new P2SoftSphericalTruncatedForceShifted(space, potentialAA, 2.5);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2TruncatedAA);
            P2LennardJones potentialAB = new P2LennardJones(space, 0.8, 1.5);
            P2SoftSphericalTruncated p2TruncatedAB = new P2SoftSphericalTruncatedForceShifted(space, potentialAB, 2.5);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), p2TruncatedAB);
            P2LennardJones potentialBB = new P2LennardJones(space, sigmaB, 0.5);
            P2SoftSphericalTruncated p2TruncatedBB = new P2SoftSphericalTruncatedForceShifted(space, potentialBB, 2.5);
            potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), p2TruncatedBB);
        } else if (potentialChoice == PotentialChoice.WCA) {
            neighborManager.setNeighborRange(2);
            // https://doi.org/10.1103/PhysRevX.1.021013
            sigmaB = D == 2 ? 1.4 : 1.0 / 1.2; // changes 1/1.4 to 1.4
            double mA = D == 2 ? 1 : 2;
            P2WCA potentialAA = new P2WCA(space, 1, 1);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), potentialAA);
            P2WCA potentialAB = new P2WCA(space, D == 2 ? 1.1 : 0.5 + 0.5 * sigmaB, 1);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), potentialAB);
            P2WCA potentialBB = new P2WCA(space, sigmaB, 1);
            potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), potentialBB);
            ((ElementSimple) speciesA.getLeafType().getElement()).setMass(mA);
        } else if (potentialChoice == PotentialChoice.SS) {
            // https://doi.org/10.1103/PhysRevLett.81.120 prescribes cut=4.5*(0.5+0.5/1.4)=3.85714
            sigmaB = 1.0 / 1.4;
            P2SoftSphere potentialAA = new P2SoftSphere(space, 1, 1, 12);
            P2SoftSphericalTruncated p2TruncatedAA = new P2SoftSphericalTruncatedForceShifted(space, potentialAA, 2.5);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2TruncatedAA);
            P2SoftSphere potentialAB = new P2SoftSphere(space, 0.5 + 0.5 * sigmaB, 1, 12);
            P2SoftSphericalTruncated p2TruncatedAB = new P2SoftSphericalTruncatedForceShifted(space, potentialAB, 2.5);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), p2TruncatedAB);
            P2SoftSphere potentialBB = new P2SoftSphere(space, sigmaB, 1, 12);
            P2SoftSphericalTruncated p2TruncatedBB = new P2SoftSphericalTruncatedForceShifted(space, potentialBB, 2.5);
            potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), p2TruncatedBB);
        } else if (potentialChoice == PotentialChoice.HS) {
            chs = 50;
            double coreHS = 0.01 * chs;
            double L = Math.pow((nA + nB) / density, 1.0 / 3.0);
            if (L < 2.01) throw new RuntimeException("too small!");
            double nbrCut = 1.7;
            if (L < nbrCut * 2) nbrCut = L / 2.001;
            neighborManager.setNeighborRange(nbrCut);
            sigmaB = 1.0 / 1.4;
            p2AA = P2SquareWell.makePotential(coreHS, 1 / coreHS, -100);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2AA);
            P2HardGeneric potentialAB = P2HardSphere.makePotential(0.5 + 0.5 * sigmaB);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), potentialAB);
            P2HardGeneric potentialBB = P2HardSphere.makePotential(sigmaB);
            potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), potentialBB);
        } else if (potentialChoice == PotentialChoice.WS) {
            // G. Wahnström, Phys. Rev. A 44, 3752 1991.
            // T. B. Schrøder, cond-mat/0005127.
            // https://doi.org/10.1063/1.1605094
            sigmaB = 5.0 / 6.0;
            P2LennardJones potentialAA = new P2LennardJones(space);
            P2SoftSphericalTruncated p2TruncatedAA = new P2SoftSphericalTruncatedForceShifted(space, potentialAA, 2.5);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2TruncatedAA);
            P2LennardJones potentialAB = new P2LennardJones(space, (1 + sigmaB) / 2, 1);
            P2SoftSphericalTruncated p2TruncatedAB = new P2SoftSphericalTruncatedForceShifted(space, potentialAB, 2.5);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), p2TruncatedAB);
            P2LennardJones potentialBB = new P2LennardJones(space, sigmaB, 1);
            P2SoftSphericalTruncated p2TruncatedBB = new P2SoftSphericalTruncatedForceShifted(space, potentialBB, 2.5);
            potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), p2TruncatedBB);
            ((ElementSimple) speciesA.getLeafType().getElement()).setMass(2);
        }

        //construct box
        box.setNMolecules(speciesA, nA);
        box.setNMolecules(speciesB, nB);
        BoxInflate boxInflate = new BoxInflate(box, box.getSpace(), density);
        boxInflate.actionPerformed();
        new ConfigurationLattice(space.D() == 2 ? (new LatticeOrthorhombicHexagonal(space)) : (new LatticeCubicFcc(space)), space).initializeCoordinates(box);

        integrator = potentialChoice == PotentialChoice.HS ?
                new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster), (NeighborManagerHard) neighborManager, random, tStep, temperature, box) :
                new IntegratorVelocityVerletFasterer(potentialMaster, random, tStep, temperature, box);
        integrator.setIsothermal(true);
        integrator.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN);
        integrator.setThermostatInterval(1);
        integrator.setThermostatNoDrift(true);
        if (potentialChoice == PotentialChoice.HS) {
            integrator.setAlwaysScaleRandomizedMomenta(true);
        }

        swapMove = new MCMoveSwapFasterer(space, random, potentialMasterMC, speciesA, speciesB);
        potentialMasterMC.setPairPotentials(potentialMaster.getPairPotentials());
        integratorMC = new IntegratorMCFasterer(potentialMasterMC, random, integrator.getTemperature(), box);
        integratorMC.getMoveManager().addMCMove(swapMove);
        if (doSwap) {
            integrator.setThermostat(IntegratorMDFasterer.ThermostatType.HYBRID_MC);
            integrator.setIntegratorMC(integratorMC, 10000);
            integrator.setThermostatInterval(1000);
        }

        integrator.reset();
        integrator.doThermostat();

        timeSource = new TimeSourceFasterer();
    }

    public Activity makeInitConfigActivity() {
        return new Activity() {
            @Override
            public void runActivity(Controller.ControllerHandle handle) {
                if (potentialChoice != PotentialChoice.HS) return;
                boolean success = false;
                PotentialComputePairGeneral potentialMaster = (PotentialComputePairGeneral) integrator.getPotentialCompute();
                PotentialComputePair potentialMasterMC = (PotentialComputePair) integratorMC.getPotentialCompute();
                double tStepOld = integrator.getTimeStep();
                integrator.setTimeStep(0.001);
                for (; chs <= 100; chs++) {
                    p2AA.setCollisionDiameter(0, chs * 0.01);
                    double u = potentialMaster.computeAll(false);
//                    System.out.println("chs "+chs*0.01+" "+u);
                    if (u == Double.POSITIVE_INFINITY) {
                        chs--;
                        p2AA.setCollisionDiameter(0, chs * 0.01);
                        break;
                    }
                    if (chs == 100) {
                        success = true;
                        P2HardGeneric p = P2HardSphere.makePotential(1);
                        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p);
                        potentialMasterMC.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p);
                    }
                }
                if (!success) {

                    while (chs < 100) {
                        integrator.reset();
                        for (int i = 0; i < 1000; i++) {
                            handle.yield(integrator::doStep);
                        }
                        chs++;
                        p2AA.setCollisionDiameter(0, chs * 0.01);
                        double u = potentialMaster.computeAll(false);
//                        System.out.println("chs "+chs*0.01+" "+u);
                        if (u == Double.POSITIVE_INFINITY) {
                            chs--;
                            p2AA.setCollisionDiameter(0, chs * 0.01);
                            continue;
                        } else if (chs == 100) {
                            P2HardGeneric p = P2HardSphere.makePotential(1);
                            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p);
                            potentialMasterMC.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p);
                        }
                    }
                }
                integrator.reset();
                integrator.resetStepCount();
                integrator.setTimeStep(tStepOld);
            }
        };
    }

    public void initConfig() {
        this.getController().runActivityBlocking(makeInitConfigActivity());
    }

    @Override
    public IntegratorMDFasterer getIntegrator() {
        return integrator;
    }

    public static void main(String[] args) {

        GlassParams params = new GlassParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
        }
        SimGlassFasterer sim = new SimGlassFasterer(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential, params.tStep);
        sim.initConfig();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));

    }//end of main

    public static class GlassParams extends ParameterBase {
        public int D = 2;
        public int nA = 130, nB = 70;
        // rho=1000/29.34^2=1.16 (yields P=0 at T=0) for LJ http://dx.doi.org/10.1088/0953-8984/21/3/035117
        //     for 3D, 1000/8.88^3 = 1.43
        //       or (with xA=0.8), rho=1.25 yields T0=1.06 (http://dx.doi.org/10.1103/PhysRevX.1.021013)
        // for SS, P=18.37 https://doi.org/10.1103/PhysRevLett.81.120
        //    rho=1.35 at T=Tg=0.55
        // for HS, rho=1.35 is very high, but 1.30 is perhaps not high enough
        // for WCA, rho=0.75*1.4*1.4 = 1.47 with T0 = 2.13 in 2D
        //          rho=1.2 with T0=0.32 in 3D
        public double density = 1000 / (29.34 * 29.34);
        public double temperature = 1.0;
        public boolean doSwap = true;
        public PotentialChoice potential = PotentialChoice.LJ;
        public double tStep = 0.005;
    }

    public class TimeSourceFasterer implements TimeSource {
        @Override
        public double getCurrentTime() {
            return integrator.getCurrentTime();
        }

        @Override
        public Box getBox() {
            return integrator.getBox();
        }

        @Override
        public long getStepCount() {
            return integrator.getStepCount();
        }

        @Override
        public double getTimeStep() {
            return integrator.getTimeStep();
        }
    }
}