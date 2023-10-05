/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.controller.Activity;
import etomica.action.controller.Controller;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.NeighborListManagerHard;
import etomica.potential.*;
import etomica.potential.compute.NeighborManagerHard;
import etomica.potential.compute.PotentialComputePair;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

public class SimGlass extends Simulation {

    public final SpeciesGeneral speciesA, speciesB;
    public final Box box;
    public final IntegratorMD integrator;

    public final MCMoveSwap swapMove;
    public final IntegratorMC integratorMC;
    public final PotentialChoice potentialChoice;
    protected P2HardGeneric p2AA, p2AB;

    public enum PotentialChoice {LJ, WCA, SS, HS, WS}

    public double sigmaB;

    protected int chs;

    public SimGlass(int D, int nA, int nB, double density, double temperature, boolean doSwap, PotentialChoice pc, double tStep, double rcLJ) {
        this(D, nA, nB, density, temperature, doSwap, pc, tStep, null, rcLJ);
    }

    public SimGlass(int D, int nA, int nB, double density, double temperature, boolean doSwap, PotentialChoice pc, double tStep, int[] randSeeds, double rcLJ) {
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

        NeighborListManager neighborManager = pc == PotentialChoice.HS
                ? new NeighborListManagerHard(getSpeciesManager(), box, 2, 2.99, BondingInfo.noBonding())
                : new NeighborListManager(getSpeciesManager(), box, 2, 1.2*rcLJ, BondingInfo.noBonding());
        NeighborCellManager neighborManagerMC = new NeighborCellManager(getSpeciesManager(), box, 2, BondingInfo.noBonding());
        PotentialComputePair potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        //construct box
        box.setNMolecules(speciesA, nA);
        box.setNMolecules(speciesB, nB);
        BoxInflate boxInflate = new BoxInflate(box, box.getSpace(), density);
        boxInflate.actionPerformed();
        new ConfigurationLattice(space.D() == 2 ? (new LatticeOrthorhombicHexagonal(space)) : (new LatticeCubicFcc(space)), space).initializeCoordinates(box);

        if (potentialChoice == PotentialChoice.LJ) { //3D KA-80-20; 2D KA-65-35
            System.out.println(" rcLJ: " + rcLJ);
            sigmaB = 0.88;
            P2LennardJones potentialAA = new P2LennardJones();
            P2SoftSphericalTruncated p2TruncatedAA = new P2SoftSphericalTruncatedForceShifted(potentialAA, rcLJ);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2TruncatedAA);
            P2LennardJones potentialAB = new P2LennardJones(0.8, 1.5);
            P2SoftSphericalTruncated p2TruncatedAB = new P2SoftSphericalTruncatedForceShifted(potentialAB, rcLJ);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), p2TruncatedAB);
            P2LennardJones potentialBB = new P2LennardJones(sigmaB, 0.5);
            P2SoftSphericalTruncated p2TruncatedBB = new P2SoftSphericalTruncatedForceShifted(potentialBB, rcLJ);
            potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), p2TruncatedBB);
        } else if (potentialChoice == PotentialChoice.WCA) {
            neighborManager.setNeighborRange(2);
            // https://doi.org/10.1103/PhysRevX.1.021013
            sigmaB = D == 2 ? 1.4 : 1.0 / 1.2; // changes 1/1.4 to 1.4
            double mA = D == 2 ? 1 : 2;
            P2WCA potentialAA = new P2WCA(1, 1);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), potentialAA);
            P2WCA potentialAB = new P2WCA(D == 2 ? 1.1 : 0.5 + 0.5 * sigmaB, 1);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), potentialAB);
            P2WCA potentialBB = new P2WCA(sigmaB, 1);
            potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), potentialBB);
            ((ElementSimple) speciesA.getLeafType().getElement()).setMass(mA);
        } else if (potentialChoice == PotentialChoice.SS) {
            // https://doi.org/10.1103/PhysRevLett.81.120 prescribes cut=4.5*(0.5+0.5/1.4)=3.85714
            sigmaB = 1.0 / 1.4;
            P2SoftSphere potentialAA = new P2SoftSphere(1, 1, 12);
            P2SoftSphericalTruncated p2TruncatedAA = new P2SoftSphericalTruncatedForceShifted(potentialAA, 2.5);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2TruncatedAA);
            P2SoftSphere potentialAB = new P2SoftSphere(0.5 + 0.5 * sigmaB, 1, 12);
            P2SoftSphericalTruncated p2TruncatedAB = new P2SoftSphericalTruncatedForceShifted(potentialAB, 2.5);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), p2TruncatedAB);
            P2SoftSphere potentialBB = new P2SoftSphere(sigmaB, 1, 12);
            P2SoftSphericalTruncated p2TruncatedBB = new P2SoftSphericalTruncatedForceShifted(potentialBB, 2.5);
            potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), p2TruncatedBB);
        } else if (potentialChoice == PotentialChoice.HS) {
            chs = getMaxCHS();
            double sigmaA = 0.01 * chs;
            double L = Math.pow((nA + nB) / density, 1.0 / 3.0);
            if (L < 2.01) throw new RuntimeException("too small!");
            double nbrCut = 1.7;
            if (L < nbrCut * 2) nbrCut = L / 2.001;
            neighborManager.setNeighborRange(nbrCut);
            sigmaB = 1.0 / 1.4;
            p2AA = P2SquareWell.makePotential(sigmaA, 1 / sigmaA, -100);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2AA);
            p2AB = P2SquareWell.makePotential(0.5*(sigmaA + sigmaB), (1+sigmaB)/(sigmaA+sigmaB), -100);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), p2AB);
            P2HardGeneric potentialBB = P2HardSphere.makePotential(sigmaB);
            potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), potentialBB);
        } else if (potentialChoice == PotentialChoice.WS) {
            // G. Wahnström, Phys. Rev. A 44, 3752 1991.
            // T. B. Schrøder, cond-mat/0005127.
            // https://doi.org/10.1063/1.1605094
            sigmaB = 5.0 / 6.0;
            P2LennardJones potentialAA = new P2LennardJones();
            P2SoftSphericalTruncated p2TruncatedAA = new P2SoftSphericalTruncatedForceShifted(potentialAA, 2.5);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2TruncatedAA);
            P2LennardJones potentialAB = new P2LennardJones((1 + sigmaB) / 2, 1);
            P2SoftSphericalTruncated p2TruncatedAB = new P2SoftSphericalTruncatedForceShifted(potentialAB, 2.5);
            potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), p2TruncatedAB);
            P2LennardJones potentialBB = new P2LennardJones(sigmaB, 1);
            P2SoftSphericalTruncated p2TruncatedBB = new P2SoftSphericalTruncatedForceShifted(potentialBB, 2.5);
            potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), p2TruncatedBB);
            ((ElementSimple) speciesA.getLeafType().getElement()).setMass(2);
        }

        integrator = potentialChoice == PotentialChoice.HS ?
                new IntegratorHard(potentialMaster.getPairPotentials(), (NeighborManagerHard) neighborManager, random, tStep, temperature, box, getSpeciesManager()) :
                new IntegratorVelocityVerlet(potentialMaster, random, tStep, temperature, box);
        integrator.setIsothermal(true);
        integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN);
        integrator.setThermostatInterval(1);
        integrator.setThermostatNoDrift(true);
        if (potentialChoice == PotentialChoice.HS) {
            integrator.setAlwaysScaleRandomizedMomenta(true);
        }

        PotentialComputePair potentialMasterMC = new PotentialComputePair(getSpeciesManager(), box, neighborManagerMC, potentialMaster.getPairPotentials());
        swapMove = new MCMoveSwap(space, random, potentialMasterMC, speciesA, speciesB);
        integratorMC = new IntegratorMC(potentialMasterMC, random, integrator.getTemperature(), box);
        integratorMC.getMoveManager().addMCMove(swapMove);
        if (doSwap) {
            integrator.setThermostat(IntegratorMD.ThermostatType.HYBRID_MC);
            integrator.setIntegratorMC(integratorMC, 10000);
            integrator.setThermostatInterval(1000);
        }

        integrator.reset();
        integrator.doThermostat();
    }

    protected int getMaxCHS() {
        double minAA = Double.POSITIVE_INFINITY, minAB = Double.POSITIVE_INFINITY;
        Vector dr = getSpace().makeVector();
        for (IAtom a1 : box.getLeafList()) {
            if (a1.getType() != speciesA.getLeafType()) continue;
            for (IAtom a2 : box.getLeafList()) {
                if (a2 == a1) continue;
                dr.Ev1Mv2(a2.getPosition(), a1.getPosition());
                box.getBoundary().nearestImage(dr);
                if (a2.getType() == speciesA.getLeafType()) minAA = Math.min(minAA, dr.squared());
                else minAB = Math.min(minAB, dr.squared());
            }
        }
        minAB = Math.sqrt(minAB);
        minAA = Math.sqrt(minAA);

        double sigmaB = 1.0/1.4;
        // 0.01 chs fSigma sigmaB < minAA
        // chs < 100 minAA / (fSigma sigmaB)
        double maxACHS = 100.0 * minAA / 1.0;

        // 0.5 (0.01 chs fSigma sigmaB + sigmaB) < minAB
        // chs < 100 (2 minAB / sigmaB - 1) / fSigma
        double maxABCHS = 100.0 * (2 * minAB / sigmaB - 1) / 1.4;
        return Math.min((int) Math.min(maxACHS, maxABCHS), 100);
    }

    public Activity makeInitConfigActivity() {
        return new Activity() {
            @Override
            public void runActivity(Controller.ControllerHandle handle) {
                if (potentialChoice != PotentialChoice.HS || chs == 100) return;
                PotentialComputePairGeneral potentialMaster = (PotentialComputePairGeneral) integrator.getPotentialCompute();
                PotentialComputePair potentialMasterMC = (PotentialComputePair) integratorMC.getPotentialCompute();
                double tStepOld = integrator.getTimeStep();
                integrator.setTimeStep(0.001);
                // run MD and increase sigma
                while (chs < 100) {
                    integrator.reset();
                    for (int i = 0; i < 1000; i++) {
                        handle.yield(integrator::doStep);
                    }
                    int newCHS = getMaxCHS();
                    if (newCHS == chs) continue;
                    chs = newCHS;
                    if (chs < 100) {
                        double sigmaA = 0.01 * chs;
                        p2AA.setCollisionDiameter(0, sigmaA);
                        p2AB.setCollisionDiameter(0, 0.5*(sigmaA + sigmaB));
                    }
                    else {
                        // potentialMaster makes a copy of integrator's potentials, so set both
                        // ironically, potentialMasterMC potentials are the same as integrator
                        P2HardGeneric pAA = P2HardSphere.makePotential(1);
                        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), pAA);
                        potentialMasterMC.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), pAA);

                        P2HardGeneric pAB = P2HardSphere.makePotential(0.5*(1 + sigmaB));
                        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), pAB);
                        potentialMasterMC.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), pAB);
                    }
                    potentialMaster.computeAll(false);

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
    public IntegratorMD getIntegrator() {
        return integrator;
    }

    public static void main(String[] args) {

        GlassParams params = new GlassParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
        }
        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential, params.tStep, params.rcLJ);
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
        public boolean doLinear = false;
        public PotentialChoice potential = PotentialChoice.LJ;
        public double tStep = 0.005;
        public double rcLJ = 2.5;
    }

}