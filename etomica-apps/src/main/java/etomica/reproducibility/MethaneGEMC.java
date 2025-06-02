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
import etomica.data.types.DataGroup;
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
        IPotential2 p2 = shift ? new P2SoftSphericalTruncatedShifted(p2lj, rc) : new P2SoftSphericalTruncated(p2lj, rc);
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
//            params.steps = 120 * 1000 * (params.numAtoms[0] + params.numAtoms[1]);
        }

        // TraPPE says CH4 is 16.04
        Unit dUnit = new SimpleUnit(Null.DIMENSION, 1/(16.04/Constants.AVOGADRO*1e24), "Density", "g/cm^3", false);
        Unit pUnit = new PrefixedUnit(Prefix.KILO, Pascal.UNIT);

        double temperatureK = params.temperatureK;
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        double[] densitygcm = params.density;
        double[] density = new double[]{dUnit.toSim(densitygcm[0]), dUnit.toSim(densitygcm[1])};

        MethaneGEMC sim = new MethaneGEMC(params.numAtoms, density, temperature, params.rc, params.shift, params.lrc);

        if (false) {
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
        long equilibration = params.equilibration;
        int interval = 10;
        int blocks = 100;
        long blockSize = Math.max(steps / (interval * blocks), 1);

        System.out.println("Lennard-Jones Monte Carlo simulation");
        System.out.println("N: " + params.numAtoms[0]+" "+params.numAtoms[1]);
        System.out.println("T: " + temperatureK);
        System.out.println("initial density: " + density[0]+" "+density[1]);
        System.out.println("initial density (g/cm^3): " + densitygcm[0]+" "+densitygcm[1]);
        System.out.println("steps: " + params.steps);

        // equilibration
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorGEMC, equilibration));
        System.out.println("equilibration finished\n");
        sim.integratorGEMC.getMoveManager().setEquilibrating(false);
        sim.integrator[0].getMoveManager().setEquilibrating(false);
        sim.integrator[1].getMoveManager().setEquilibrating(false);
        sim.mcMoveAtom[0].getTracker().reset();
        sim.mcMoveAtom[1].getTracker().reset();

        // data collection
        MeterPotentialEnergyFromIntegrator meterPEliq = new MeterPotentialEnergyFromIntegrator(sim.integrator[0]);
        AccumulatorAverageFixed accPEliq = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpPEliq = new DataPumpListener(meterPEliq, accPEliq, interval);
        sim.integrator[0].getEventManager().addListener(pumpPEliq);

        MeterPotentialEnergyFromIntegrator meterPEvap = new MeterPotentialEnergyFromIntegrator(sim.integrator[1]);
        AccumulatorAverageFixed accPEvap = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpPEvap = new DataPumpListener(meterPEvap, accPEvap, interval);
        sim.integrator[1].getEventManager().addListener(pumpPEvap);

        MeterDensity meterDensityLiq = new MeterDensity(sim.integrator[0].getBox());
        AccumulatorAverageFixed accDensityLiq = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpDensityLiq = new DataPumpListener(meterDensityLiq, accDensityLiq, interval);
        sim.integrator[0].getEventManager().addListener(pumpDensityLiq);

        MeterDensity meterDensityVap = new MeterDensity(sim.integrator[1].getBox());
        AccumulatorAverageFixed accDensityVap = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpDensityVap = new DataPumpListener(meterDensityVap, accDensityVap, interval);
        sim.integrator[1].getEventManager().addListener(pumpDensityVap);

        sim.integratorGEMC.resetStepCount();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorGEMC, steps));

        long t2 = System.currentTimeMillis();

        DataGroup dataPEliq = (DataGroup) accPEliq.getData();
        int numAtomsLiq = sim.getBox(0).getLeafList().size();
        double avgLiq = dataPEliq.getValue(accPEliq.AVERAGE.index) / numAtomsLiq;
        double errLiq = dataPEliq.getValue(accPEliq.ERROR.index) / numAtomsLiq;
        double corLiq = dataPEliq.getValue(accPEliq.BLOCK_CORRELATION.index);

        DataGroup dataPEvap = (DataGroup) accPEvap.getData();
        int numAtomsVap = sim.getBox(0).getLeafList().size();
        double avgVap = dataPEvap.getValue(accPEvap.AVERAGE.index) / numAtomsVap;
        double errVap = dataPEvap.getValue(accPEvap.ERROR.index) / numAtomsVap;
        double corVap = dataPEvap.getValue(accPEvap.BLOCK_CORRELATION.index);

        DataGroup dataDensityLiq = (DataGroup) accDensityLiq.getData();
        double avgDensityLiq = dataDensityLiq.getValue(accDensityLiq.AVERAGE.index);
        double errDensityLiq = dataDensityLiq.getValue(accDensityLiq.ERROR.index);
        double corDensityLiq = dataDensityLiq.getValue(accDensityLiq.BLOCK_CORRELATION.index);

        DataGroup dataDensityVap = (DataGroup) accDensityVap.getData();
        double avgDensityVap = dataDensityVap.getValue(accDensityVap.AVERAGE.index);
        double errDensityVap = dataDensityVap.getValue(accDensityVap.ERROR.index);
        double corDensityVap = dataDensityVap.getValue(accDensityVap.BLOCK_CORRELATION.index);

        System.out.println("liquid energy avg: " + avgLiq + "  err: " + errLiq + "  cor: " + corLiq);
        System.out.println("vapor energy avg: " + avgVap + "  err: " + errVap + "  cor: " + corVap);
        System.out.println("liquid density avg: " + avgDensityLiq + "  err: " + errDensityLiq + "  cor: " + corDensityLiq);
        System.out.println("vapor density avg: " + avgDensityVap + "  err: " + errDensityVap + "  cor: " + corDensityVap);
        System.out.println("liquid density avg (g/cm^3): " + dUnit.fromSim(avgDensityLiq) + "  err: " + dUnit.fromSim(errDensityLiq) + "  cor: " + corDensityLiq);
        System.out.println("vapor density avg (g/cm^3): " + dUnit.fromSim(avgDensityVap) + "  err: " + dUnit.fromSim(errDensityVap) + "  cor: " + corDensityVap);
        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class SimParams extends ParameterBase {
        public int[] numAtoms = new int[]{900,100};
        public long steps = (numAtoms[0]+numAtoms[1])*1000*120;
        public long equilibration = (numAtoms[0]+numAtoms[1])*1000*50;
        public double[] density = new double[]{0.3752,0.0117};
        public double temperatureK = 140;
        public double rc = 14;
        public boolean shift = false;
        public boolean lrc = true;
    }
}
