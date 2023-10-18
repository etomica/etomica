/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.controller.Activity;
import etomica.action.controller.Controller;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterPressureHardA;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.NeighborListManagerHard;
import etomica.potential.BondingInfo;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2HardSphere;
import etomica.potential.P2SquareWell;
import etomica.potential.compute.PotentialComputePair;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import javax.swing.*;
import java.awt.event.WindowEvent;

/**
 * This simulation measures the free energy derivative for growing the diameter
 * of a single atom in a system of atoms that are equal or smaller sizes.
 */
public class SimSolidHSA1 extends Simulation {


    public final NeighborListManagerHard neighborManager;
    public final SpeciesGeneral speciesA, speciesB;
    public final Box box;
    public final IntegratorHard integrator;
    protected P2HardGeneric p2AA, p2AB;
    public double sigmaB, fSigma, chs;

    public SimSolidHSA1(int D, int n, double fSigma, double density, double temperature, double tStep) {
        this(D, n, fSigma, density, temperature, tStep, null);
    }

    public SimSolidHSA1(int D, int n, double fSigma, double density, double temperature, double tStep, int[] randSeeds) {
        super(Space.getInstance(D));
        if (randSeeds != null) {
            setRandom(new RandomMersenneTwister(randSeeds));
        }

        this.fSigma = fSigma;

        //species
        speciesA = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesA);
        speciesB = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesB);

        box = this.makeBox();

        neighborManager = new NeighborListManagerHard(getSpeciesManager(), box, 2, 2.99, BondingInfo.noBonding());
        PotentialComputePair potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        double L = Math.pow(n / density, 1.0 / 3.0);
        if (L < 2.01) throw new RuntimeException("too small!");
        double nbrCut = 1.7;
        if (L < nbrCut * 2) nbrCut = L / 2.001;
        neighborManager.setNeighborRange(nbrCut);
        sigmaB = 1.0 / 1.4;

        double rnbr = Math.pow(4/density, 1.0/3.0)/Math.sqrt(2);
        chs = 100;
        double sigmaA = 0.01 * chs * fSigma * sigmaB;
        double sigmaAB = 0.5*(sigmaA + sigmaB);
        while (sigmaAB > rnbr) {
            chs--;
            sigmaA = 0.01 * chs * fSigma * sigmaB;
            sigmaAB = 0.5*(sigmaA + sigmaB);
        }

        p2AA = P2SquareWell.makePotential(sigmaA, fSigma * sigmaB / sigmaA, -100);
        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2AA);
        p2AB = P2SquareWell.makePotential(0.5*(sigmaA + sigmaB), (fSigma*sigmaB + sigmaB)/(sigmaA+sigmaB), -100);
        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), p2AB);

        P2HardGeneric potentialBB = P2HardSphere.makePotential(sigmaB);
        potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), potentialBB);

        //construct box
        box.setNMolecules(speciesA, 1);
        box.setNMolecules(speciesB, n-1);
        BoxInflate boxInflate = new BoxInflate(box, box.getSpace(), density);
        boxInflate.actionPerformed();
        new ConfigurationLattice(space.D() == 2 ? (new LatticeOrthorhombicHexagonal(space)) : (new LatticeCubicFcc(space)), space).initializeCoordinates(box);

        integrator = new IntegratorHard(potentialMaster.getPairPotentials(), neighborManager, random, tStep, temperature, box, getSpeciesManager());
        integrator.setIsothermal(false);
        integrator.setThermostatNoDrift(true);
        integrator.setAlwaysScaleRandomizedMomenta(true);

        integrator.reset();
        integrator.doThermostat();
    }

    public Activity makeInitConfigActivity() {
        return new Activity() {
            @Override
            public void runActivity(Controller.ControllerHandle handle) {
                if (chs == 100) return;
                PotentialComputePairGeneral potentialMaster = (PotentialComputePairGeneral) integrator.getPotentialCompute();
                double tStepOld = integrator.getTimeStep();
                integrator.setTimeStep(0.001);
                // if sigma still too small, run MD and continue increasing sigma
                while (chs < 100) {
                    integrator.reset();
                    for (int i = 0; i < 1000; i++) {
                        handle.yield(integrator::doStep);
                    }
                    chs++;
                    double sigmaA = 0.01 * chs * fSigma * sigmaB;
                    double sigmaAB = 0.5*(sigmaA + sigmaB);
                    p2AA.setCollisionDiameter(0, sigmaA);
                    p2AB.setCollisionDiameter(0, sigmaAB);
                    double u = potentialMaster.computeAll(false);
                    if (u == Double.POSITIVE_INFINITY) {
                        chs--;
                        sigmaA = 0.01 * chs * fSigma * sigmaB;
                        sigmaAB = 0.5*(sigmaA + sigmaB);
                        p2AA.setCollisionDiameter(0, sigmaA);
                        p2AA.setCollisionDiameter(1, 1/sigmaA);
                        p2AB.setCollisionDiameter(0, sigmaAB);
                        p2AB.setCollisionDiameter(1, 0.5*(fSigma + 1));
                    } else if (chs == 100) {
                        P2HardGeneric pAA = P2HardSphere.makePotential(1);
                        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), pAA);

                        P2HardGeneric pAB = P2HardSphere.makePotential(0.5*(1 + sigmaB));
                        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), pAB);
                    }
                }
                integrator.reset();
                integrator.resetStepCount();
                integrator.setTimeStep(tStepOld);
            }
        };
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
            params.D = 3;
            params.n = 256;
            params.density = 2.9;
            params.fSigma = 1.01;
            params.steps = 1000000;
        }

        System.out.println(params.D+"D HS B, growing single A");
        System.out.println("N: "+params.n);
        System.out.println("fSigmaA: "+params.fSigma);
        System.out.println("density: "+params.density);
        System.out.println();

        SimSolidHSA1 sim = new SimSolidHSA1(params.D, params.n, params.fSigma, params.density, params.temperature, params.tStep, params.seeds);

        if (false) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "HS Fluid A1");
            DiameterHashByType diameterHash = (DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash();
            diameterHash.setDiameter(sim.speciesA.getLeafType(), params.fSigma*sim.sigmaB);
            diameterHash.setDiameter(sim.speciesB.getLeafType(), sim.sigmaB);
            simGraphic.makeAndDisplayFrame();

            JFrame f = new JFrame();
            f.setSize(700, 500);
            JPanel panel = new JPanel();
            panel.add(new JLabel("<html><div style='width: 200px;'>Note: high-density simulations will be slightly delayed in starting up, as the simulation works to generate a configuration without overlaps. Controls probably won't work properly until this dialog disappears. </div></html>"));
            f.getContentPane().add(panel);
            f.pack();
            f.setTitle("Generating configuration");
            f.setLocationRelativeTo(null);
            f.setVisible(true);
            sim.getController().addActivity(sim.makeInitConfigActivity()).future.whenComplete((res, ex) -> {
                simGraphic.getController().getResetAveragesButton().getAction().actionPerformed();
                f.dispatchEvent(new WindowEvent(f, WindowEvent.WINDOW_CLOSING));
            });

            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));

            return;
        }

        sim.getController().runActivityBlocking(sim.makeInitConfigActivity());

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.steps/10));
        System.out.println("Equilibration finished");
        sim.integrator.resetStepCount();

        int interval = 10;
        long bs = params.steps/(interval * 100);
        if (bs == 0) bs = 1;

        MeterPressureHardA meterPA = new MeterPressureHardA(sim.integrator, sim.speciesA.getLeafType());
        AccumulatorAverageFixed accPA = new AccumulatorAverageFixed(bs);
        DataPumpListener pumpPA = new DataPumpListener(meterPA, accPA, interval);
        sim.integrator.getEventManager().addListener(pumpPA);

        MeterPressureHard meterP = new MeterPressureHard(sim.integrator);
        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(bs);
        DataPumpListener pumpP = new DataPumpListener(meterP, accP, interval);
        sim.integrator.getEventManager().addListener(pumpP);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.steps));

        IData dataPA = accPA.getData();
        double avgPA = dataPA.getValue(AccumulatorAverageFixed.AVERAGE.index);
        double errPA = dataPA.getValue(AccumulatorAverageFixed.ERROR.index);
        double corPA = dataPA.getValue(AccumulatorAverageFixed.BLOCK_CORRELATION.index);
        System.out.println("PA: "+avgPA+"  err: "+errPA+"  cor: "+corPA);

        // dFdsigmaA = 0.5 * PA * 3 * V / sigmaA
        // 0.5 is because sigmaAB increases half as fast as sigmaA

        IData dataP = accP.getData();
        double avgP = dataP.getValue(AccumulatorAverageFixed.AVERAGE.index);
        double errP = dataP.getValue(AccumulatorAverageFixed.ERROR.index);
        double corP = dataP.getValue(AccumulatorAverageFixed.BLOCK_CORRELATION.index);
        System.out.println("P: "+avgP+"  err: "+errP+"  cor: "+corP);

        long t2 = System.nanoTime();
        System.out.println("time: "+(t2-t1)/1e9);

    }//end of main

    public static class GlassParams extends ParameterBase {
        public int D = 3;
        public int n = 200;
        public double density = 1.2;
        public double temperature = 1.0;
        public double tStep = 0.01;

        public long steps = 1000000;
        public int[] seeds;
        public double fSigma = 1;
    }

}