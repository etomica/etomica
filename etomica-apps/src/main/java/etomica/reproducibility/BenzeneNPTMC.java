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
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataDouble;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveMoleculeRotate;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBenzene;
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
public class BenzeneNPTMC extends Simulation {

    public PotentialCompute pcAgg;
    public IntegratorMC integrator;
    public SpeciesGeneral species;
    public Box box;
    public MCMoveVolume mcMoveVolume;
    public MCMoveMolecule translateMove;

    public BenzeneNPTMC(Space space, double density, int nSpheres, int numMolecules, double temperature, double rc, double pressure) {
        super(space);
        setRandom(new RandomMersenneTwister(1));
        species = SpeciesBenzene.makeBuilder().build();
        addSpecies(species);

        box = this.makeBox();
        box.setNMolecules(species, numMolecules);
        new BoxInflate(box, space, density).actionPerformed();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);

        double k = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT).toSim(938);
        P2Harmonic p2Bond = new P2Harmonic(Kelvin.UNIT.toSim(k), SpeciesBenzene.nominalBondL);
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<nSpheres; i++) {
            int[] p = new int[]{i,i+1};
            for (int j=0; j<p.length; j++) {
                if (p[j] >= 6) p[j] -= 6;
            }
            pairs.add(p);
        }
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        P3BondAngle p3Bond = new P3BondAngle(Math.PI*120.0/180.0, Kelvin.UNIT.toSim(70450));
        List<int[]> triplets = new ArrayList<>();
        for (int i=0; i<nSpheres; i++) {
            int[] t = new int[]{i,i+1,i+2};
            for (int j=0; j<t.length; j++) {
                if (t[j] >= 6) t[j] -= 6;
            }
            triplets.add(t);
        }
        pmBonding.setBondingPotentialTriplet(species, p3Bond, triplets);

        P4BondTorsionPow p4Bond = new P4BondTorsionPow(new double[]{Kelvin.UNIT.toSim(260), 0, Kelvin.UNIT.toSim(-530), 0,
                0,0,0,0, Kelvin.UNIT.toSim(530)});
        List<int[]> quads = new ArrayList<>();
        for (int i=0; i<nSpheres; i++) {
            int[] q = new int[]{i,i+1,i+2,i+3};
            for (int j=0; j<q.length; j++) {
                if (q[j] >= 6) q[j] -= 6;
            }
            quads.add(q);
        }
        pmBonding.setBondingPotentialQuad(species, p4Bond, quads);

        PotentialMasterCell potentialMaster = new PotentialMasterCell(getSpeciesManager(), box, 2, pmBonding.getBondingInfo());
        potentialMaster.doAllTruncationCorrection = true;
        pcAgg = new PotentialComputeAggregate(pmBonding, potentialMaster);
        integrator = new IntegratorMC(pcAgg, random, temperature, box);

        translateMove = new MCMoveMolecule(random, pcAgg, box);
        integrator.getMoveManager().addMCMove(translateMove);

        MCMoveMoleculeRotate rotateMove = new MCMoveMoleculeRotate(random, pcAgg, box);
        integrator.getMoveManager().addMCMove(rotateMove);

        MCMoveAtom moveAtom = new MCMoveAtom(random, pcAgg, box);
        integrator.getMoveManager().addMCMove(moveAtom);


        if (pressure >= 0) {
            mcMoveVolume = new MCMoveVolume(integrator, random, pressure);
            integrator.getMoveManager().addMCMove(mcMoveVolume);
        }

        AtomType typeCH = species.getAtomType(0);
        double epsilon = Kelvin.UNIT.toSim(50.5);
        double sigma = 3.695;
        TruncationFactory tf = new TruncationFactorySimple(rc);
        P2LennardJones p2LJ = new P2LennardJones(sigma, epsilon);
        IPotential2 p2 = tf.make(p2LJ);

        potentialMaster.setPairPotential(typeCH, typeCH, p2);


        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
        potentialMaster.init();
        double u0 = potentialMaster.computeAll(false);
        double x = 1;
        while (u0 > 1e7*numMolecules) {
            x *= 0.99;
            p2LJ.setSigma(x*sigma);
            ((P2SoftSphericalSumTruncatedForceShifted)p2).setTruncationRadius(rc);
            u0 = potentialMaster.computeAll(false);
        }
        integrator.reset();
        while (u0 > 1e7*numMolecules) {
            double d = box.getMoleculeList().size() / box.getBoundary().volume();
            System.out.println(integrator.getStepCount()+" "+u0/numMolecules+" "+x+" "+d);
            while (u0 > 1e7 * numMolecules) {
                integrator.doStep();
                u0 = integrator.getPotentialEnergy();
                System.out.println(integrator.getStepCount()+" "+u0/numMolecules+" "+x+" "+d);
            }
            while (x < 1 && u0 <= 1e7*numMolecules) {
                x /= 0.99;
                if (x > 1) x = 1;
                p2LJ.setSigma(x*sigma);
                ((P2SoftSphericalSumTruncatedForceShifted)p2).setTruncationRadius(rc);
                u0 = potentialMaster.computeAll(false);
            }
            integrator.reset();
        }
    }

    public static void main(String[] args) {

        BenzeneParams params = new BenzeneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.numSteps = 350000;
            params.density = 0.3;
            params.graphics = true;
        }

        Unit dUnit = new SimpleUnit(Null.DIMENSION, 1/(78.11/Constants.AVOGADRO*1e24), "Density", "g/cm^3", false);

        double temperatureK = params.temperatureK;
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        System.out.println("Tsim "+temperature);
        int numMolecules = params.numMolecules;
        double density = dUnit.toSim(params.density);
        boolean graphics = params.graphics;
        long numSteps = params.numSteps;
        double rc = params.rc;
        double pressureKPa = params.pressureKPa;
        Unit pUnit = new PrefixedUnit(Prefix.KILO, Pascal.UNIT);
        double pressure = pUnit.toSim(pressureKPa);

        System.out.println(numSteps+" steps");
        System.out.println("rc: "+rc);
        System.out.println("pressure "+ pressureKPa);
        System.out.println("initial density "+ density);
        System.out.println("initial density (g/cm^3) "+ dUnit.fromSim(density));

        final BenzeneNPTMC sim = new BenzeneNPTMC(Space3D.getInstance(), density, 5, numMolecules, temperature, rc, pressure);

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

            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Benzene MC", 3);

            DiameterHashByType dhbt = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
            dhbt.setDiameter(sim.species.getAtomType(0), 3.695);

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

            simGraphic.makeAndDisplayFrame();
            return;
        }

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps/100));
        sim.integrator.getMoveManager().addMCMove(sim.translateMove);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps/5));

        long samples = numSteps / (numMolecules*8);
        long bs = samples / 100;
        if (bs == 0) bs = 1;

        AccumulatorAverageFixed accU = new AccumulatorAverageFixed((numSteps/10)/100);
        DataPumpListener pumpU = new DataPumpListener(meterU, accU, 10);
        sim.integrator.getEventManager().addListener(pumpU);

        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(bs);
        forkP.addDataSink(accP);
        AccumulatorAverageFixed accZ = new AccumulatorAverageFixed(bs);
        dpZ.addDataSink(accZ);
        AccumulatorAverageFixed accZm1oR = new AccumulatorAverageFixed(bs);
        dpZm1oR.addDataSink(accZm1oR);
        DataPumpListener pumpP = new DataPumpListener(meterP, forkP, 8*numMolecules);
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

        System.out.println("time: "+(t2-t1)/1e9);
    }

    public static class BenzeneParams extends ParameterBase {
        public double temperatureK = 450;
        public int numMolecules = 400;
        public double density = 0.692;
        public boolean graphics = false;
        public long numSteps = 200000;
        public double rc = 14;
        public double pressureKPa = 2260;
    }
}
