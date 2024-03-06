/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.reproducibility;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveMoleculeRotate;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.simulation.prototypes.MCMoveWiggle;
import etomica.simulation.prototypes.MeterTorsionAngle;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesAlkane;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.util.ArrayList;
import java.util.List;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class PentaneNPTMC extends Simulation {

    public PotentialCompute pcAgg;
    public IntegratorMC integrator;
    public SpeciesGeneral species;
    public Box box;
    public MCMoveVolume mcMoveVolume;
    public MCMoveMolecule translateMove;

    public PentaneNPTMC(Space space, double density, int nSpheres, int numMolecules, double temperature, String configFilename, double rc, double pressure) {
        super(space);
        setRandom(new RandomMersenneTwister(1));
        species = SpeciesAlkane.makeBuilder(nSpheres)
                .build();
        addSpecies(species);

        box = this.makeBox();
        box.setNMolecules(species, numMolecules);
        new BoxInflate(box, space, density).actionPerformed();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);

        P3BondAngle p3 = new P3BondAngle(Math.PI*114.0/180.0, Kelvin.UNIT.toSim(62500));
        List<int[]> triplets = new ArrayList<>();
        for (int i=0; i<nSpheres-2; i++) {
            triplets.add(new int[]{i,i+1,i+2});
        }
        pmBonding.setBondingPotentialTriplet(species, p3, triplets);

        P4BondTorsion p4 = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
        List<int[]> quads = new ArrayList<>();
        for (int i=0; i<nSpheres-3; i++) {
            quads.add(new int[]{i,i+1,i+2,i+3});
        }
        pmBonding.setBondingPotentialQuad(species, p4, quads);

        PotentialMasterCell potentialMaster = new PotentialMasterCell(getSpeciesManager(), box, 2, pmBonding.getBondingInfo());
        potentialMaster.doAllTruncationCorrection = true;
        pcAgg = new PotentialComputeAggregate(pmBonding, potentialMaster);
        integrator = new IntegratorMC(pcAgg, random, temperature, box);

        translateMove = new MCMoveMolecule(random, pcAgg, box);
        integrator.getMoveManager().addMCMove(translateMove);

        MCMoveMoleculeRotate rotateMove = new MCMoveMoleculeRotate(random, pcAgg, box);
        integrator.getMoveManager().addMCMove(rotateMove);

        MCMoveWiggle wiggleMove = new MCMoveWiggle(random, pcAgg, box);
        integrator.getMoveManager().addMCMove(wiggleMove);

        if (pressure >= 0) {
            mcMoveVolume = new MCMoveVolume(integrator, random, pressure);
            integrator.getMoveManager().addMCMove(mcMoveVolume);
        }

        AtomType typeCH3 = species.getAtomType(0);
        AtomType typeCH2 = species.getAtomType(1);
        double epsilonCH2 = Kelvin.UNIT.toSim(46.0);
        double epsilonCH3 = Kelvin.UNIT.toSim(98.0);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
        double sigmaCH2 = 3.95;
        double sigmaCH3 = 3.75;
        double sigmaCH2CH3 = (sigmaCH2+sigmaCH3)/2;
        TruncationFactory tf = new TruncationFactoryForceShift(rc);
        P2LennardJones p2CH2LJ = new P2LennardJones(sigmaCH2, epsilonCH2);
        P2LennardJones p2CH3LJ = new P2LennardJones(sigmaCH3, epsilonCH3);
        P2LennardJones p2CH2CH3LJ = new P2LennardJones(sigmaCH2CH3, epsilonCH2CH3);
        IPotential2 p2CH2 = tf.make(p2CH2LJ);
        IPotential2 p2CH3 = tf.make(p2CH3LJ);
        IPotential2 p2CH2CH3 = tf.make(p2CH2CH3LJ);

        potentialMaster.setPairPotential(typeCH2, typeCH2, p2CH2);
        potentialMaster.setPairPotential(typeCH2, typeCH3, p2CH2CH3);
        potentialMaster.setPairPotential(typeCH3, typeCH3, p2CH3);

        if (configFilename != null) {
            ConfigurationFile config = new ConfigurationFile(configFilename);
            config.initializeCoordinates(box);
            BoxImposePbc.imposePBC(box);
        }
        else {
            ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
            configuration.initializeCoordinates(box);
            potentialMaster.init();
            double u0 = potentialMaster.computeAll(false);
            double x = 1;
            while (u0 > 1e6*numMolecules) {
                x *= 0.99;
                p2CH2LJ.setSigma(x*sigmaCH2);
                p2CH3LJ.setSigma(x*sigmaCH3);
                p2CH2CH3LJ.setSigma(x*sigmaCH2CH3);
                ((P2SoftSphericalSumTruncatedForceShifted)p2CH2).setTruncationRadius(rc);
                ((P2SoftSphericalSumTruncatedForceShifted)p2CH3).setTruncationRadius(rc);
                ((P2SoftSphericalSumTruncatedForceShifted)p2CH2CH3).setTruncationRadius(rc);
                u0 = potentialMaster.computeAll(false);
            }
            integrator.reset();
            while (u0 > 1e4*numMolecules) {
                while (u0 > 1e4 * numMolecules) {
                    integrator.doStep();
                    u0 = integrator.getPotentialEnergy();
                }
                while (x < 1 && u0 <= 1e4*numMolecules) {
                    x /= 0.99;
                    if (x > 1) x = 1;
                    p2CH2LJ.setSigma(x*sigmaCH2);
                    p2CH3LJ.setSigma(x*sigmaCH3);
                    p2CH2CH3LJ.setSigma(x*sigmaCH2CH3);
                    ((P2SoftSphericalSumTruncatedForceShifted)p2CH2).setTruncationRadius(rc);
                    ((P2SoftSphericalSumTruncatedForceShifted)p2CH3).setTruncationRadius(rc);
                    ((P2SoftSphericalSumTruncatedForceShifted)p2CH2CH3).setTruncationRadius(rc);
                    u0 = potentialMaster.computeAll(false);
                }
                integrator.reset();
            }
        }
    }

    public static void main(String[] args) {

        OctaneParams params = new OctaneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.numSteps = 350000;
            params.density = 0.3;
            params.configFilename = null; // "octane";
            params.graphics = false;
        }

        Unit dUnit = new SimpleUnit(Null.DIMENSION, 1/(72.15/Constants.AVOGADRO*1e24), "Density", "g/cm^3", false);

        double temperatureK = params.temperatureK;
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        System.out.println("Tsim "+temperature);
        int numMolecules = params.numMolecules;
        double density = dUnit.toSim(params.density);
        boolean graphics = params.graphics;
        long numSteps = params.numSteps;
        String configFilename = params.configFilename;
        double rc = params.rc;
        double pressureKPa = params.pressureKPa;
        Unit pUnit = new PrefixedUnit(Prefix.KILO, Pascal.UNIT);
        double pressure = pUnit.toSim(pressureKPa);

        System.out.println(numSteps+" steps");
        System.out.println("rc: "+rc);
        System.out.println("pressure "+ pressureKPa);
        System.out.println("initial density "+ density);
        System.out.println("initial density (g/cm^3) "+ dUnit.fromSim(density));

        final PentaneNPTMC sim = new PentaneNPTMC(Space3D.getInstance(), density, 5, numMolecules, temperature, configFilename, rc, pressure);

        MeterPotentialEnergyFromIntegrator meterU = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        sim.integrator.getPotentialCompute().init();
        sim.integrator.reset();
        System.out.println("u0/N "+(meterU.getDataAsScalar()/numMolecules));

        MeterPressure meterP = new MeterPressure(sim.box, sim.pcAgg);
        meterP.setTemperature(temperature);
        meterP.doCallComputeAll(true);
        DataProcessorForked dpZ = new DataProcessorForked() {
            DataDouble.DataInfoDouble dataInfo = new DataDouble.DataInfoDouble("Z", Null.DIMENSION);
            DataDouble data = new DataDouble();

            @Override
            protected IData processData(IData inputData) {
                data.x = inputData.getValue(0) / temperature / density;
                return data;
            }

            @Override
            protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                return dataInfo;
            }
        };
        DataProcessorForked dpZm1oR = new DataProcessorForked() {
            DataDouble.DataInfoDouble dataInfo = new DataDouble.DataInfoDouble("(Z-1)/rho", Null.DIMENSION);
            DataDouble data = new DataDouble();

            @Override
            protected IData processData(IData inputData) {
                data.x = (inputData.getValue(0) / temperature / density - 1) / density;
                return data;
            }

            @Override
            protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                return dataInfo;
            }
        };
        DataFork forkP = new DataFork(new IDataSink[]{dpZ, dpZm1oR});

        if (graphics) {
//            sim.getController().addActivity(new ActivityIntegrate(sim.integrator, numMolecules*100));
//            sim.getController().addActionSequential(new IAction() {
//                @Override
//                public void actionPerformed() {
//                    sim.integrator.getMoveManager().addMCMove(sim.mcMoveVolume);
//                }
//            });
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));

            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Pentane MC", 3);

            DiameterHashByType dhbt = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
            dhbt.setDiameter(sim.species.getAtomType(0), 3.75);
            dhbt.setDiameter(sim.species.getAtomType(1), 3.95);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            List<DataPump> dataPumps = simGraphic.getController().getDataStreamPumps();

            DataSourceCountSteps timer = new DataSourceCountSteps(sim.integrator);
            DisplayTextBox timerBox = new DisplayTextBox();
            timerBox.setLabel("Steps");
            DataPumpListener pumpSteps = new DataPumpListener(timer, timerBox, numMolecules);
            sim.integrator.getEventManager().addListener(pumpSteps);
            simGraphic.add(timerBox);

            MeterDensity meterDensity = new MeterDensity(sim.box());
            AccumulatorHistory accDensity = new AccumulatorHistory(new HistoryCollapsingAverage());
            accDensity.setTimeDataSource(timer);
            DataPumpListener pumpDensity = new DataPumpListener(meterDensity, accDensity, 10);
            sim.integrator.getEventManager().addListener(pumpDensity);
            dataPumps.add(pumpDensity);

            DisplayPlot historyDensity = new DisplayPlot();
            accDensity.setDataSink(historyDensity.getDataSet().makeDataSink());
            historyDensity.setLabel("Density");
            historyDensity.setUnit(dUnit);
            simGraphic.add(historyDensity);

            Unit perN = new SimpleUnit(Null.DIMENSION, numMolecules, "1/N", "1/N", false);

            AccumulatorHistory historyU = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyU.setTimeDataSource(timer);
            AccumulatorHistory historyU2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyU2.setTimeDataSource(timer);
            AccumulatorAverageCollapsing avgEnergy = new AccumulatorAverageCollapsing();
            avgEnergy.setPushInterval(10);
            DataFork forkU = new DataFork(new IDataSink[]{historyU, historyU2, avgEnergy});
            DataPumpListener pumpU = new DataPumpListener(meterU, forkU, 10);
            dataPumps.add(pumpU);
            sim.integrator.getEventManager().addListener(pumpU);
            DisplayPlotXChart plotU = new DisplayPlotXChart();
            plotU.setLabel("U");
            historyU.addDataSink(plotU.makeSink("U"));
            plotU.setLegend(new DataTag[]{historyU.getTag()}, "samples");
            historyU2.addDataSink(plotU.makeSink("Uavg"));
            plotU.setLegend(new DataTag[]{historyU2.getTag()}, "avg");
            plotU.setUnit(perN);
            simGraphic.add(plotU);

            simGraphic.getController().getDataStreamPumps().add(pumpU);

            DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
            display.setAccumulator(avgEnergy);
            display.setUnit(perN);
            simGraphic.add(display);

            meterP.setTemperature(temperature);
            meterP.doCallComputeAll(true);
            AccumulatorHistory historyP = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyP.setTimeDataSource(timer);
            AccumulatorHistory historyP2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyP2.setTimeDataSource(timer);
            AccumulatorAverageCollapsing avgP = new AccumulatorAverageCollapsing();
            forkP.addDataSink(historyP);
            forkP.addDataSink(historyP2);
            forkP.addDataSink(avgP);
            DataPumpListener pumpP = new DataPumpListener(meterP, forkP, numMolecules);
            dataPumps.add(pumpP);
            sim.integrator.getEventManager().addListener(pumpP);
            DisplayPlotXChart plotP = new DisplayPlotXChart();
            plotP.setLabel("P");
            historyP.addDataSink(plotP.makeSink("P"));
            plotP.setLegend(new DataTag[]{historyP.getTag()}, "samples");
            historyP2.addDataSink(plotP.makeSink("Pavg"));
            plotP.setLegend(new DataTag[]{historyP2.getTag()}, "avg");
            plotP.setUnit(pUnit);
            simGraphic.add(plotP);
            simGraphic.getController().getDataStreamPumps().add(pumpP);

            DisplayTextBoxesCAE displayP = new DisplayTextBoxesCAE();
            displayP.setAccumulator(avgP);
            simGraphic.add(displayP);

            AccumulatorHistory historyZ = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyZ.setTimeDataSource(timer);
            AccumulatorHistory historyZ2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyZ2.setTimeDataSource(timer);
            AccumulatorAverageCollapsing avgZ = new AccumulatorAverageCollapsing();
            dpZ.addDataSink(historyZ);
            dpZ.addDataSink(historyZ2);
            dpZ.addDataSink(avgZ);
            DisplayPlotXChart plotZ = new DisplayPlotXChart();
            plotZ.setLabel("Z");
            historyZ.addDataSink(plotZ.makeSink("Zsamples"));
            plotZ.setLegend(new DataTag[]{historyZ.getTag()}, "samples");
            historyZ2.addDataSink(plotZ.makeSink("Zavg"));
            plotZ.setLegend(new DataTag[]{historyZ2.getTag()}, "avg");
            simGraphic.add(plotZ);

            DisplayTextBoxesCAE displayZ = new DisplayTextBoxesCAE();
            displayZ.setLabel("Z");
            displayZ.setAccumulator(avgZ);
            simGraphic.add(displayZ);

            AccumulatorHistory historyZ_ = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyZ_.setTimeDataSource(timer);
            AccumulatorHistory historyZ_2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyZ_2.setTimeDataSource(timer);
            AccumulatorAverageCollapsing avgZ_ = new AccumulatorAverageCollapsing();
            dpZm1oR.addDataSink(historyZ_);
            dpZm1oR.addDataSink(historyZ_2);
            dpZm1oR.addDataSink(avgZ_);
            DisplayPlotXChart plotZ_ = new DisplayPlotXChart();
            plotZ_.setLabel("(Z-1)/rho");
            historyZ_.addDataSink(plotZ_.makeSink("samples"));
            plotZ_.setLegend(new DataTag[]{historyZ_.getTag()}, "samples");
            historyZ_2.addDataSink(plotZ_.makeSink("avg"));
            plotZ_.setLegend(new DataTag[]{historyZ_2.getTag()}, "avg");
            simGraphic.add(plotZ_);

            DisplayTextBoxesCAE displayZ_ = new DisplayTextBoxesCAE();
            displayZ_.setLabel("(Z-1)/rho");
            displayZ_.setAccumulator(avgZ_);
            simGraphic.add(displayZ_);

            MeterTorsionAngle meterTorsion = new MeterTorsionAngle(sim.box, 2, 3, 4, 5);
            AccumulatorAverageBlockless accTorsion1 = new AccumulatorAverageBlockless();
            AccumulatorAverageCollapsing accTorsion2 = new AccumulatorAverageCollapsing();
            AccumulatorHistory historyTorsion = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyTorsion.setTimeDataSource(timer);
            AccumulatorHistory historyTorsion2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyTorsion2.setTimeDataSource(timer);
            DataFork forkTorsion = new DataFork(new IDataSink[]{accTorsion1, accTorsion2, historyTorsion, historyTorsion2});
            DataPumpListener pumpTorsion = new DataPumpListener(meterTorsion, forkTorsion, numMolecules);
            dataPumps.add(pumpTorsion);
            sim.integrator.getEventManager().addListener(pumpTorsion);
            DisplayTextBoxesCAE displayTorsion = new DisplayTextBoxesCAE();
            displayTorsion.setLabel("Torsion cos");
            displayTorsion.setAccumulator(accTorsion2);
            DisplayPlotXChart plotTorsion = new DisplayPlotXChart();
            plotTorsion.setLabel("torsion");
            historyTorsion.addDataSink(plotTorsion.makeSink("samples"));
            plotTorsion.setLegend(new DataTag[]{historyTorsion.getTag()}, "samples");
            historyTorsion2.addDataSink(plotTorsion.makeSink("avg"));
            plotTorsion.setLegend(new DataTag[]{historyTorsion2.getTag()}, "avg");
            simGraphic.add(plotTorsion);
            simGraphic.makeAndDisplayFrame();
            return;
        }

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps/5));

        int interval = 10;
        int pInterval = 10*numMolecules;
        long samples = numSteps / interval;
        long bs = Math.max(samples / 100, 1);
        long pSamples = numSteps / pInterval;
        long pbs = Math.max(pSamples / 100, 1);

        MeterDensity meterDensity = new MeterDensity(sim.box());
        AccumulatorAverageFixed accDensity = new AccumulatorAverageFixed(bs);
        DataPumpListener pumpDensity = new DataPumpListener(meterDensity, accDensity, interval);
        sim.integrator.getEventManager().addListener(pumpDensity);

        AccumulatorAverageFixed accU = new AccumulatorAverageFixed(bs);
        DataPumpListener pumpU = new DataPumpListener(meterU, accU, interval);
        sim.integrator.getEventManager().addListener(pumpU);

        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(pbs);
        forkP.addDataSink(accP);
        AccumulatorAverageFixed accZ = new AccumulatorAverageFixed(pbs);
        dpZ.addDataSink(accZ);
        AccumulatorAverageFixed accZm1oR = new AccumulatorAverageFixed(pbs);
        dpZm1oR.addDataSink(accZm1oR);
        DataPumpListener pumpP = new DataPumpListener(meterP, forkP, pInterval);
        sim.integrator.getEventManager().addListener(pumpP);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));
        long t2 = System.nanoTime();

        IData dataU = accU.getData();
        double avgU = dataU.getValue(accU.AVERAGE.index) / numMolecules;
        double errU = dataU.getValue(accU.ERROR.index) / numMolecules;
        double corU = dataU.getValue(accU.BLOCK_CORRELATION.index);
        System.out.println("U: "+avgU+"   err: "+errU+"   cor: "+corU);

        IData dataP = accP.getData();
        double avgP = dataP.getValue(accP.AVERAGE.index);
        double errP = dataP.getValue(accP.ERROR.index);
        double corP = dataP.getValue(accP.BLOCK_CORRELATION.index);
        System.out.println("P: "+avgP+"   err: "+errP+"   cor: "+corP);

        IData dataZ = accZ.getData();
        double avgZ = dataZ.getValue(accZ.AVERAGE.index);
        double errZ = dataZ.getValue(accZ.ERROR.index);
        double corZ = dataZ.getValue(accZ.BLOCK_CORRELATION.index);
        System.out.println("Z: "+avgZ+"   err: "+errZ+"   cor: "+corZ);

        IData dataZ_ = accZm1oR.getData();
        double avgZ_ = dataZ_.getValue(accZm1oR.AVERAGE.index);
        double errZ_ = dataZ_.getValue(accZm1oR.ERROR.index);
        double corZ_ = dataZ_.getValue(accZm1oR.BLOCK_CORRELATION.index);
        System.out.println("(Z-1)/rho: "+avgZ_+"   err: "+errZ_+"   cor: "+corZ_);

        DataGroup dataDensity = (DataGroup) accDensity.getData();
        double avgDensity = dataDensity.getValue(accDensity.AVERAGE.index);
        double errDensity = dataDensity.getValue(accDensity.ERROR.index);
        double corDensity = dataDensity.getValue(accDensity.BLOCK_CORRELATION.index);
        System.out.println("density avg: " + avgDensity + "  err: " + errDensity + "  cor: " + corDensity);
        System.out.println("density avg (g/cm^3): " + dUnit.fromSim(avgDensity) + "  err: " + dUnit.fromSim(errDensity) + "  cor: " + corDensity);

        System.out.println("time: "+(t2-t1)/1e9);
    }

    public static class OctaneParams extends ParameterBase {
        public double temperatureK = 372;
        public int numMolecules = 300;
        public double density = 0.0005;
        public boolean graphics = false;
        public long numSteps = 200000;
        public String configFilename = null;
        public double rc = 14;
        public double pressureKPa = 1402;
    }
}
