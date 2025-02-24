/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicBcc;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.BondingInfo;
import etomica.potential.P1Sinusoidal;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Simple Lennard-Jones Monte Carlo simulation in 3D.
 */
public class SimLJLiquidTI extends Simulation {
    public IntegratorMC integrator;
//    public IntegratorLangevin integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesGeneral species;
    public Box box;
    public PotentialMasterCell pmMeter;
    public PotentialComputeField pc1Meter;
    public P2LennardJones potential;
    public double temperature;

    public SimLJLiquidTI(int nAtoms, double rho, double temperature, double rc, double kSine, double lambda, Structure strc, boolean isLRC) {
        super(Space3D.getInstance());
//        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        AtomType leafType = species.getLeafType();
        addSpecies(species);
        box = this.makeBox();
        box.setNMolecules(species, nAtoms);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(rho);
        inflater.actionPerformed();

        ConfigurationLattice config;
        String strcString;
        if (strc == Structure.SC) {
            config = new ConfigurationLattice(new LatticeCubicSimple(space), space);
            strcString = "SC";
        } else if (strc == Structure.BCC) {
            config = new ConfigurationLattice(new LatticeCubicBcc(space), space);
            strcString = "BCC";
        } else { //FCC
            config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
            strcString = "FCC";
        }
        config.initializeCoordinates(box);

        double sigma = 1.0;
        double eps = 1.0;
        System.out.println(" L/2: " + box.getBoundary().getBoxSize().getX(0)/2.0);

        potential = new P2LennardJones(sigma, eps*lambda);
        P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(potential, rc);
        PotentialMasterCell  pm = new PotentialMasterCell(getSpeciesManager(), box, 2, BondingInfo.noBonding());
        pm.setPairPotential(leafType, leafType, potentialTruncated);

        P2LennardJones potentialMeter = new P2LennardJones(sigma, eps);
        P2SoftSphericalTruncated potentialTruncatedMeter = new P2SoftSphericalTruncated(potentialMeter, rc);
        pmMeter = new PotentialMasterCell(getSpeciesManager(), box, 2, BondingInfo.noBonding());
        pmMeter.setPairPotential(leafType, leafType, potentialTruncatedMeter);

        if (!isLRC) {
            pm.doAllTruncationCorrection = false;
            pmMeter.doAllTruncationCorrection = false;
        }

        double a;
        if (strc == Structure.SC) {
            a = Math.pow(1.0/rho, 1.0/3.0);
        } else if (strc == Structure.BCC) {
            a = Math.pow(2.0/rho, 1.0/3.0);
        } else {
            a = Math.pow(4.0/rho, 1.0/3.0);
        }
        System.out.println(" a: " + a);

        Vector shift = space.makeVector();
        shift.E(box.getLeafList().get(0).getPosition());

        P1Sinusoidal p1Sinusoidal = new P1Sinusoidal(getSpace(), a, (1-lambda)*kSine, strcString, shift);
        PotentialComputeField pc1 = new PotentialComputeField(getSpeciesManager(), box);
        pc1.setFieldPotential(leafType, p1Sinusoidal);

        PotentialComputeAggregate.localStorageDefault = true;
        PotentialComputeAggregate pcAgg = new PotentialComputeAggregate(pc1, pm);
//        double gamma = 10;
//        integrator = new IntegratorLangevin(pcAgg, random, 0.001, temperature, box, gamma);
//        integrator.setThermostatNoDrift(false);
//        integrator.setIsothermal(true);

        integrator = new IntegratorMC(pcAgg, random, temperature, box);
        mcMoveAtom = new MCMoveAtom(random, pcAgg, box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        P1Sinusoidal p1SinusoidalMeter = new P1Sinusoidal(getSpace(), a, kSine, strcString, shift);
        pc1Meter = new PotentialComputeField(getSpeciesManager(), box);
        pc1Meter.setFieldPotential(leafType, p1SinusoidalMeter);

        // init() are required. Note needed for pm and pc1 because the integrator calls it.
//        pmMeter.init();
//        pc1Meter.init();
    }

    public static void main(String[] args) {
        long t1 = System.nanoTime();
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int nAtoms = params.nAtoms;
        int nSteps = params.nSteps;
        int nStepsEq = params.nSteps/10;
        double temperature = params.temperature;
        double rho = params.rho;
        double rc = params.rc;
        Structure strc = params.strc;
        double lambda = params.lambda;
        double kSine = params.kSine;
        boolean isGraphic = params.isGraphic;
        boolean isLRC = params.isLRC;

        SimLJLiquidTI sim = new SimLJLiquidTI(nAtoms, rho, temperature, rc, kSine, lambda, strc, isLRC);
        System.out.println(" LJ");
        System.out.println(" N: " + nAtoms);
        System.out.println(" density: " + rho);
        System.out.println(" T: " + temperature);
        System.out.println(" rc: " + rc);
        System.out.println(" isLRC: " + isLRC);
        System.out.println(" steps: " +  nSteps);
        System.out.println(" kSine: " + kSine);
        System.out.println(" lambda: " + lambda);
        System.out.println(" Structure: " + strc);

        if (isGraphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.setPaintInterval(sim.box, 1000);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

                public Color getAtomColor(IAtom a) {
                    if (allColors==null) {
                        allColors = new Color[768];
                        for (int i=0; i<256; i++) {
                            allColors[i] = new Color(255-i,i,0);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+256] = new Color(0,255-i,i);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+512] = new Color(i,0,255-i);
                        }
                    }
                    return allColors[(2*a.getLeafIndex()) % 768];
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

            ((DiameterHashByType)((DisplayBox)simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species().getAtomType(0),1);

            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());

            DisplayPlot pePlot = new DisplayPlot();
            simGraphic.add(pePlot);

            simGraphic.makeAndDisplayFrame();
            return;
        }
        System.out.flush();

        // Equilibaration
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, nStepsEq));
        System.out.println(" Done with " + nStepsEq + " steps equilibration.");

        int interval = nAtoms;
        int blocks = 100;
        long blockSize = params.nSteps / (interval * blocks);
        System.out.println(" nBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);
        System.out.println();

        sim.pmMeter.init();
        sim.pc1Meter.init();
//        MeterPotentialEnergyFromIntegrator meterDUdlambda = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        MeterdUdlambda meterDUdlambda = new MeterdUdlambda(sim.pmMeter, sim.pc1Meter);
        AccumulatorAverageFixed accumulatorDUdlambda = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpDUdlambda = new DataPumpListener(meterDUdlambda, accumulatorDUdlambda, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpDUdlambda);


        // Production
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, nSteps));
        System.out.println(" Move acceptance: " + sim.mcMoveAtom.getTracker().acceptanceProbability());

        double avg = accumulatorDUdlambda.getData(accumulatorDUdlambda.AVERAGE).getValue(0);
        double err = accumulatorDUdlambda.getData(accumulatorDUdlambda.ERROR).getValue(0);
        double cor = accumulatorDUdlambda.getData(accumulatorDUdlambda.BLOCK_CORRELATION).getValue(0);
        System.out.println();
        System.out.println(" dUdl: " + avg/nAtoms + "  " + err/nAtoms + "    " + cor);

        long t2 = System.nanoTime();
        System.out.println("\n time: " + (t2 - t1)/1.0e9/60.0 + " min");

    }

    public enum Structure {SC, BCC, FCC};

    public static class SimParams extends ParameterBase {
        public Structure strc = Structure.SC;
        public boolean isGraphic = false;
        public boolean isLRC = !false;
        int n = 7;
        public int nAtoms = 1*n*n*n;
        public int nSteps = 1000000;
        public double temperature = 1;
        public double rho = 0.8;
        public double rc = 3;
        public double kSine = 0;
        public double lambda = 1;
    }
}