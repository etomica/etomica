/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.reproducibility;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.*;
import etomica.graphics.DisplayPlotXChart;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorGEMC;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.function.Function;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.*;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Simulation for Monte Carlo simulation of Lennard-Jonesium.  Interactions are
 * tracked with cell lists.
 *
 * @author Andrew Schultz
 */

public class MethaneGEMC extends Simulation {

    public final PotentialMasterCell[] potentialMaster;
    public final IntegratorManagerMC integratorGEMC;
    public final IntegratorMC[] integrator;
    public MCMoveAtom[] mcMoveAtom;

    /**
     * Creates simulation with the given parameters
     */
    public MethaneGEMC(int[] numAtoms, double[] density, double temperature, double rc, boolean shift, boolean lrc) {
        super(Space3D.getInstance());

        SpeciesGeneral species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        double epsilon = Kelvin.UNIT.toSim(148);
        double sigma = 3.73;
        P2LennardJones p2lj = new P2LennardJones(sigma, epsilon);
        Potential2Soft p2 = shift ? new P2SoftSphericalTruncatedShifted(p2lj, rc) : new P2SoftSphericalTruncated(p2lj, rc);
        AtomType atomType = species.getLeafType();

        potentialMaster = new PotentialMasterCell[2];
        integrator = new IntegratorMC[2];
        mcMoveAtom = new MCMoveAtom[2];
        Box[] boxes = new Box[2];
        for (int i=0; i<2; i++) {
            boxes[i] = new Box(space);
            addBox(boxes[i]);
            potentialMaster[i] = new PotentialMasterCell(getSpeciesManager(), boxes[i], 2, BondingInfo.noBonding());
            potentialMaster[i].doAllTruncationCorrection = lrc;
            potentialMaster[i].doOneTruncationCorrection = lrc;

            boxes[i].setNMolecules(species, numAtoms[i]);
            new BoxInflate(boxes[i], space, density[i]).actionPerformed();

            new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(boxes[i]);

            potentialMaster[i].setPairPotential(atomType, atomType, p2);

            integrator[i] = new IntegratorMC(potentialMaster[i], random, temperature, boxes[i]);

            mcMoveAtom[i] = new MCMoveAtom(random, potentialMaster[i], boxes[i]);
            integrator[i].getMoveManager().addMCMove(mcMoveAtom[i]);
        }

        integratorGEMC = IntegratorGEMC.buildGEMC(integrator[0], integrator[1], random, space);

        getController().addActivity(new ActivityIntegrate(integratorGEMC));
    }

    public static void main(String[] args) {

        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // modify parameters here for interactive testing
            params.steps = 120 * 1000 * (params.numAtoms[0] + params.numAtoms[1]);
        }

        // TraPPE says CH4 is 16.04
        Unit dUnit = new SimpleUnit(Null.DIMENSION, 1/(16.04/Constants.AVOGADRO*1e24), "Density", "g/cm^3", false);
        Unit pUnit = new PrefixedUnit(Prefix.KILO, Pascal.UNIT);

        double temperatureK = params.temperatureK;
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        double[] density = new double[]{dUnit.toSim(params.density[0]), dUnit.toSim(params.density[1])};

        MethaneGEMC sim = new MethaneGEMC(params.numAtoms, density, temperature, params.rc, params.shift, params.lrc);

        if (true) {
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DiameterHashByType diameterHash = (DiameterHashByType)graphic.getDisplayBox(sim.getBox(0)).getDiameterHash();
            diameterHash.setDiameter(sim.species().getLeafType(), 3.73);

            String[] l = new String[]{"L","V"};

            DisplayPlotXChart historyPressure = new DisplayPlotXChart();
            historyPressure.setLabel("pressure");
            historyPressure.setUnit(pUnit);
            DisplayPlotXChart historyVolume = new DisplayPlotXChart();
            historyVolume.setLabel("volume");
            DisplayPlotXChart historyNMolecules = new DisplayPlotXChart();
            historyNMolecules.setLabel("#");
            int N = params.numAtoms[0] + params.numAtoms[1];

            for (int i=0; i<2; i++) {
                DataSourceCountSteps timeSource = new DataSourceCountSteps(sim.integrator[i]);
                MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator[i]);
                int finalI = i;
                DataProcessorFunction pe1oN = new DataProcessorFunction(new Function() {
                    @Override
                    public double f(double x) {
                        return x/sim.integrator[finalI].getBox().getLeafList().size();
                    }
                });
                DataPumpListener pumpPE = new DataPumpListener(meterPE, pe1oN, 10);
                AccumulatorHistory accPE = new AccumulatorHistory(new HistoryCollapsingAverage());
                accPE.setTimeDataSource(timeSource);
                pe1oN.setDataSink(accPE);
                sim.integrator[i].getEventManager().addListener(pumpPE);

                DisplayPlotXChart historyPE = new DisplayPlotXChart();
                accPE.setDataSink(historyPE.getDataSet().makeDataSink());
                historyPE.setLabel("PE " + l[i]);
                graphic.add(historyPE);

                MeterDensity meterDensity = new MeterDensity(sim.getBox(i));
                AccumulatorHistory accDensity = new AccumulatorHistory(new HistoryCollapsingAverage());
                accDensity.setTimeDataSource(timeSource);
                DataPumpListener pumpDensity = new DataPumpListener(meterDensity, accDensity, 10);
                sim.integrator[i].getEventManager().addListener(pumpDensity);

                DisplayPlotXChart historyDensity = new DisplayPlotXChart();
                accDensity.setDataSink(historyDensity.getDataSet().makeDataSink());
                historyDensity.setLabel("Density");
                historyDensity.setUnit(dUnit);
                graphic.add(historyDensity);

                MeterVolume meterVolume = new MeterVolume();
                meterVolume.setBox(sim.getBox(i));
                AccumulatorHistory accVolume = new AccumulatorHistory(new HistoryCollapsingAverage());
                accVolume.setTimeDataSource(timeSource);
                DataPumpListener pumpVolume = new DataPumpListener(meterVolume, accVolume, 10);
                sim.integrator[i].getEventManager().addListener(pumpVolume);

                accVolume.setDataSink(historyVolume.getDataSet().makeDataSink());
                historyVolume.setLegend(new DataTag[]{accVolume.getTag()}, l[i]);

                MeterNMolecules meterNMolecules = new MeterNMolecules();
                meterNMolecules.setBox(sim.getBox(i));
                AccumulatorHistory accNMolecules = new AccumulatorHistory(new HistoryCollapsingAverage());
                accNMolecules.setTimeDataSource(timeSource);
                DataPumpListener pumpNMolecules = new DataPumpListener(meterNMolecules, accNMolecules, 10);
                sim.integrator[i].getEventManager().addListener(pumpNMolecules);

                accNMolecules.setDataSink(historyNMolecules.getDataSet().makeDataSink());
                historyNMolecules.setLegend(new DataTag[]{accNMolecules.getTag()}, l[i]);

                MeterPressure meterPressure = new MeterPressure(sim.getBox(i), sim.potentialMaster[i]);
                meterPressure.setTemperature(temperature);
                AccumulatorHistory accPressure = new AccumulatorHistory(new HistoryCollapsingAverage());
                accPressure.setTimeDataSource(timeSource);
                DataPumpListener pumpPressure = new DataPumpListener(meterPressure, accPressure, N);
                sim.integrator[i].getEventManager().addListener(pumpPressure);

                accPressure.setDataSink(historyPressure.getDataSet().makeDataSink());
                historyPressure.setLegend(new DataTag[]{accPressure.getTag()}, l[i]);
            }
            graphic.add(historyPressure);
            graphic.add(historyNMolecules);
            graphic.add(historyVolume);

            graphic.makeAndDisplayFrame();
            return;
        }

        long steps = params.steps;
        int interval = 10;
        int blocks = 100;
        long blockSize = Math.max(steps / (interval * blocks), 1);

        System.out.println("Lennard-Jones Monte Carlo simulation");
        System.out.println("N: " + params.numAtoms[0]+" "+params.numAtoms[1]);
        System.out.println("T: " + temperatureK);
        System.out.println("initial density: " + params.density[0]+" "+params.density[1]);
        System.out.println("steps: " + params.steps);

        // equilibration
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorGEMC, steps / 10));
        System.out.println("equilibration finished\n");
        sim.integratorGEMC.getMoveManager().setEquilibrating(false);
        sim.integrator[0].getMoveManager().setEquilibrating(false);
        sim.integrator[1].getMoveManager().setEquilibrating(false);
        sim.mcMoveAtom[0].getTracker().reset();
        sim.mcMoveAtom[1].getTracker().reset();

        // data collection
//        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
//        AccumulatorAverageFixed accPE = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE, interval);
//        sim.integrator.getEventManager().addListener(pumpPE);
//
//        MeterDensity meterDensity = new MeterDensity(sim.box());
//        AccumulatorAverageFixed accDensity = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener pumpDensity = new DataPumpListener(meterDensity, accDensity, interval);
//        sim.integrator.getEventManager().addListener(pumpDensity);

        sim.integratorGEMC.resetStepCount();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorGEMC, steps));

        long t2 = System.currentTimeMillis();
//
//        MCMoveStepTracker volumeTracker = (MCMoveStepTracker) sim.mcMoveVolume.getTracker();
//        System.out.println("volume change fraction: "+((double)volumeTracker.nTrials)/params.steps);
//        System.out.println("volume change step size: "+sim.mcMoveVolume.getStepSize());
//        System.out.println("volume change acceptance: "+volumeTracker.acceptanceProbability());
//
//        MCMoveStepTracker displacementTracker = (MCMoveStepTracker) sim.mcMoveAtom.getTracker();
//        System.out.println("displacement fraction: "+((double)displacementTracker.nTrials)/params.steps);
//        System.out.println("displacement step size: "+sim.mcMoveAtom.getStepSize());
//        System.out.println("displacement acceptance: "+displacementTracker.acceptanceProbability());
//        System.out.println();
//
//        DataGroup dataPE = (DataGroup) accPE.getData();
//        int numAtoms = sim.getBox(0).getLeafList().size();
//        double avg = dataPE.getValue(accPE.AVERAGE.index) / numAtoms;
//        double err = dataPE.getValue(accPE.ERROR.index) / numAtoms;
//        double cor = dataPE.getValue(accPE.BLOCK_CORRELATION.index);
//
//        DataGroup dataDensity = (DataGroup) accDensity.getData();
//        double avgDensity = dataDensity.getValue(accPE.AVERAGE.index);
//        double errDensity = dataDensity.getValue(accPE.ERROR.index);
//        double corDensity = dataDensity.getValue(accPE.BLOCK_CORRELATION.index);
//
//        System.out.println("energy avg: " + avg + "  err: " + err + "  cor: " + cor);
//        System.out.println("density avg: " + avgDensity + "  err: " + errDensity + "  cor: " + corDensity);
//        System.out.println("density avg (g/cm^3): " + dUnit.fromSim(avgDensity) + "  err: " + dUnit.fromSim(errDensity) + "  cor: " + corDensity);
        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class SimParams extends ParameterBase {
        public long steps = 100000000;
        public double[] density = new double[]{0.3752,0.0117};
        public double temperatureK = 140;
        public int[] numAtoms = new int[]{900,100};
        public double rc = 14;
        public boolean shift = false;
        public boolean lrc = true;
    }
}
