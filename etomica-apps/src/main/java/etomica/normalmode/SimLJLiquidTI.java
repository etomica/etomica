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
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicBcc;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.list.NeighborListManager;
import etomica.potential.BondingInfo;
import etomica.potential.P1Sinusoidal;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.numerical.BesselFunction;

import java.awt.*;

/**
 * Simple Lennard-Jones Monte Carlo simulation in 3D.
 */
public class SimLJLiquidTI extends Simulation {
    public IntegratorMC integrator;
//    public IntegratorVelocityVerlet integrator;

    public MCMoveAtom mcMoveAtom;
    public SpeciesGeneral species;
    public Box box;
    public PotentialComputePair pmMeter;
    public PotentialComputeField pc1Meter;
    public P2LennardJones potential;
    public double temperature;

    public SimLJLiquidTI(int nAtoms, double rho, double temperature, double rc, double kSine, double lambda, Structure strc) {
        super(Space3D.getInstance());
//        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);
        box = this.makeBox();
        int cellRange = 2;
        NeighborCellManager neighborManager = new NeighborCellManager(getSpeciesManager(), box, cellRange, BondingInfo.noBonding());
//        NeighborListManager neighborManager = new NeighborListManager(getSpeciesManager(), box, 2, 4, BondingInfo.noBonding());

        PotentialComputePair pm = new PotentialComputePair(getSpeciesManager(), box, neighborManager);
        pmMeter = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

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
        System.out.println(" L/2: " + box.getBoundary().getBoxSize().getX(0)/2.0);

        double sigma = 1.0;
        double eps = 1.0;
        potential = new P2LennardJones(sigma, eps*lambda);
        P2SoftSphericalTruncatedShifted potentialTruncated = new P2SoftSphericalTruncatedShifted(potential, rc);
        AtomType leafType = species.getLeafType();
        pm.setPairPotential(leafType, leafType, potentialTruncated);
        pm.doAllTruncationCorrection = false;

        P2LennardJones potentialMeter = new P2LennardJones(1,1);
        P2SoftSphericalTruncatedShifted potentialTruncatedMeter = new P2SoftSphericalTruncatedShifted(potentialMeter, rc);
        pmMeter.setPairPotential(leafType, leafType, potentialTruncatedMeter);
        pmMeter.doAllTruncationCorrection = false;

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

//        integrator = new IntegratorVelocityVerlet(pcAgg, random, 0.001, temperature, box);
//        integrator.setIsothermal(true);
//        integrator.setThermostatNoDrift(true);
//
        integrator = new IntegratorMC(pcAgg, random, temperature, box);
        mcMoveAtom = new MCMoveAtom(random, pcAgg, box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        P1Sinusoidal p1SinusoidalMeter = new P1Sinusoidal(getSpace(), a, kSine, strcString, shift);
        pc1Meter = new PotentialComputeField(getSpeciesManager(), box);
        pc1Meter.setFieldPotential(species.getLeafType(), p1SinusoidalMeter);

        double x = kSine/temperature/2;
        double Aref = -getSpace().D()*Math.log(a*Math.exp(-x)*BesselFunction.I(0, x));
        System.out.println(" Aref: " + Aref);
    }

    public static void main(String[] args) {
        long t1 = System.nanoTime();
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int nAtoms = params.nAtoms;
        int nSteps = params.nSteps;
        int nStepsEq = nSteps/10;
        double temperature = params.temperature;
        double rho = params.rho;
        double rc = params.rc;
        boolean isGraphic = params.isGraphic;
        Structure strc = params.strc;
        double lambda = params.lambda;
        double kSine = params.kSine;

        SimLJLiquidTI sim = new SimLJLiquidTI(nAtoms, rho, temperature, rc, kSine, lambda, strc);
        System.out.println(" LJ");
        System.out.println(" N: " + nAtoms);
        System.out.println(" density: " + rho);
        System.out.println(" T: " + temperature);
        System.out.println(" rc: " + rc);
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

            ((DiameterHashByType)((DisplayBox)simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species().getAtomType(0),0.7);

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
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.nSteps / 10));
        System.out.println(" Done with " + nStepsEq + " steps equilibration.");

        int interval = nAtoms;
        int blocks = 100;
        long blockSize = params.nSteps / (interval * blocks);
        System.out.println(" nBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);
        System.out.println();

        MeterdUdlambda meterDUdlambda = new MeterdUdlambda(sim.pmMeter, sim.pc1Meter);
        AccumulatorAverageFixed accumulatorDUdlambda = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpDUdlambda = new DataPumpListener(meterDUdlambda, accumulatorDUdlambda, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpDUdlambda);

        // Production
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.nSteps));
        System.out.println(" Move acceptance: " + sim.mcMoveAtom.getTracker().acceptanceProbability());

        double avg = accumulatorDUdlambda.getData(accumulatorDUdlambda.AVERAGE).getValue(0);
        double err = accumulatorDUdlambda.getData(accumulatorDUdlambda.ERROR).getValue(0);
        double cor = accumulatorDUdlambda.getData(accumulatorDUdlambda.BLOCK_CORRELATION).getValue(0);
        System.out.println();
        System.out.println(" dUdl: " + avg/nAtoms + "  " + err/nAtoms + "    " + cor);
        System.out.println( "3/2 kT: " + 3.0/2.0*temperature);


        long t2 = System.nanoTime();
        System.out.println("\n time: " + (t2 - t1)/1.0e9/60.0 + " min");
    }

    public enum Structure {SC, BCC, FCC};

    public static class SimParams extends ParameterBase {
        public Structure strc = Structure.SC;
        public boolean isGraphic = true;
        int n = 7;
        public int nAtoms = 1*n*n*n;
        public int nSteps = 1000000;
        public double temperature = 1;
        public double rho = 0.8;
        public double rc = 2.5;
        public double kSine = 5;
        public double lambda = 0;
    }
}