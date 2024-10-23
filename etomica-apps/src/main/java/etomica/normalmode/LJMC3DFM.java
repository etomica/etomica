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
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.BondingInfo;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.NeighborManagerSimple;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;

/**
 * Simple Lennard-Jones Monte Carlo simulation in 3D.
 */
public class LJMC3DFM extends Simulation {

    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesGeneral species;
    public Box box;
    public PotentialComputePair potentialMaster;
    //    public PotentialMasterList potentialMaster;
    public P2LennardJones potential;
    public double temperature;

    public LJMC3DFM(int numAtoms, double rho, double temperature, double rc) {
        super(Space3D.getInstance());








        setRandom(new RandomMersenneTwister(10));




















        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);
        box = this.makeBox();
        int cellRange = 2;
        NeighborCellManager neighborManager = new NeighborCellManager(getSpeciesManager(), box, cellRange, BondingInfo.noBonding());
        potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        double sigma = 1.0;
        double eps = 1.0;
        integrator = new IntegratorMC(potentialMaster, random, temperature, box);
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        box.setNMolecules(species, numAtoms);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(rho);
        inflater.actionPerformed();
        potential = new P2LennardJones(sigma, eps);
        P2SoftSphericalTruncatedForceShifted potentialTruncated = new P2SoftSphericalTruncatedForceShifted(potential, rc);
        AtomType leafType = species.getLeafType();
        potentialMaster.setPairPotential(leafType, leafType, potentialTruncated);

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        long t1 = System.nanoTime();
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        int numSteps = params.numSteps;
        double alpha = params.alpha;
        double temperature = params.temperature;
        double rho = params.rho;
        double rc = params.rc;

        LJMC3DFM sim = new LJMC3DFM(numAtoms, rho, temperature, rc);
        System.out.println(" LJ");
        System.out.println(" N: " + numAtoms);
        System.out.println(" density: " + rho);
        System.out.println(" T: " + temperature);
        System.out.println(" rc: " + rc);
        System.out.println(" steps: " +  numSteps);
        System.out.println(" alpha: " + alpha);

        boolean isGraphic = !true;
        if (isGraphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY);
            simGraphic.setPaintInterval(sim.box, 1);

            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

                public Color getAtomColor(IAtom a) {
                    if (allColors == null) {
                        allColors = new Color[768];
                        for (int i = 0; i < 256; i++) {
                            allColors[i] = new Color(255 - i, i, 0);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 256] = new Color(0, 255 - i, i);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 512] = new Color(i, 0, 255 - i);
                        }
                    }
                    return allColors[(768 * a.getIndex())];
                }
            };

            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species().getAtomType(0), 1);
            simGraphic.makeAndDisplayFrame("LJMCFM");

            return;
        }
        System.out.flush();

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / 10));

        int interval = numAtoms;
        int blocks = 100;
        long blockSize = params.numSteps / (interval * blocks);
        System.out.println(" numBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);

        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(blockSize);
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator, interval);
        sim.integrator.getEventManager().addListener(energyPump);

        MeterFM energyFM = new MeterFM(sim.potentialMaster, temperature, sim.box, alpha);
        AccumulatorAverage energyFMAccumulator = new AccumulatorAverageFixed(blockSize);

        DataPumpListener energyFMPump = new DataPumpListener(energyFM, energyFMAccumulator, interval);
        sim.integrator.getEventManager().addListener(energyFMPump);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));
        System.out.println("Move acceptance: " + sim.mcMoveAtom.getTracker().acceptanceProbability());

//        double avg = energyAccumulator.getData(energyAccumulator.AVERAGE).getValue(0)/numAtoms;
//        double err = energyAccumulator.getData(energyAccumulator.ERROR).getValue(0)/numAtoms;
//        double cor = energyAccumulator.getData(energyAccumulator.BLOCK_CORRELATION).getValue(0);
//        System.out.println("e_int: " + avg + " " + err + " " + cor);


//        double avg_conv = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(0);
//        double err_conv = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(0);
//        double cor_conv = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(0);
//        System.out.println("e_conv: " + avg_conv + " " + err_conv + " " + cor_conv);
//
//        double avg_fm = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(1);
//        double err_fm = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(1);
//        double cor_fm = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(1);
//        System.out.println("e_fm: " + avg_fm + " " + err_fm + " " + cor_fm);

//        double avg_fm2 = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(2);
//        double err_fm2 = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(2);
//        double cor_fm2 = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(2);
//        System.out.println("diff: " + avg_fm2 + " " + err_fm2 + " " + cor_fm2);


        double sumF2 = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(0);
        double sumLambda = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(1);
        double sumLambda2 = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(2);
        double sumFHF = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(3);
        System.out.println("sumF2  +  sumLambda  +  sumLambda2  +  sumFHF");
        System.out.println(sumF2 +" "+ sumLambda +" "+ sumLambda2 +" "+ sumFHF);

        double beta1 = 1.0/temperature + 0.01;
        double alphaOpt = (beta1*sumF2-sumLambda)/(beta1*sumFHF+sumLambda2);
        double beta0 = 1.0/temperature;
        double dAlphaOpt = sumF2/(beta0*sumFHF+sumLambda2);
        System.out.println("alphaOpt: " + alphaOpt);
        System.out.println("dAlphaOpt: " + dAlphaOpt);

        double avg_conv = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(4);
        double err_conv = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(4);
        double cor_conv = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(4);
        System.out.println("<du>: " + avg_conv + " " + err_conv + " " + cor_conv);

        double avg_fm = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(5);
        double err_fm = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(5);
        double cor_fm = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(5);
        System.out.println("<phi>: " + avg_fm + " " + err_fm + " " + cor_fm);

        double avg_dphi = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(6);
        double err_dphi = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(6);
        double cor_dphi = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(6);
        System.out.println("d<phi>/dalpha: " + avg_dphi + " " + err_dphi + " " + cor_dphi);

        double avg = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(7);
        double err = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(7);
        double cor = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(7);
        System.out.println("U0: " + avg + " " + err + " " + cor);

        long t2 = System.nanoTime();

        System.out.println("\ntime: " + (t2 - t1)/1.0e9/60.0 + " min");
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 256;
        public int numSteps = 200000;
        public double temperature = 0.4;
        public double rho = 0.9;
        public double rc = 2.5;
        public double alpha = 1e-5;
    }
}