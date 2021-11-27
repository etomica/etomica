/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.reproducibility;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByElement;
import etomica.box.Box;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationMinimized;
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
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class EthanolNPTMC extends Simulation {

    public PotentialCompute pcAgg;
    public IntegratorMC integrator;
    public SpeciesGeneral species;
    public Box box;
    public MCMoveVolume mcMoveVolume;
    public MCMoveMolecule translateMove;

    public EthanolNPTMC(Space space, double density, int nSpheres, int numMolecules, double temperature, double rc, double pressure) {
        super(space);

        List<ConformationMinimized.BondInfo> bonds = new ArrayList<>();
        bonds.add(new ConformationMinimized.BondInfo(0, 1, 1.529));
        bonds.add(new ConformationMinimized.BondInfo(1, 2, 1.41));
        bonds.add(new ConformationMinimized.BondInfo(0, 3, 1.09));
        bonds.add(new ConformationMinimized.BondInfo(0, 4, 1.09));
        bonds.add(new ConformationMinimized.BondInfo(0, 5, 1.09));
        bonds.add(new ConformationMinimized.BondInfo(1, 6, 1.09));
        bonds.add(new ConformationMinimized.BondInfo(1, 7, 1.09));
        bonds.add(new ConformationMinimized.BondInfo(2, 8, 0.945));
        ConformationMinimized conf = new ConformationMinimized(bonds, random);

        AtomType typeCT135 = AtomType.element(Carbon.INSTANCE, "CT135");
        AtomType typeCT157 = AtomType.element(Carbon.INSTANCE, "CT157");
        AtomType typeOH = AtomType.element(Oxygen.INSTANCE, "OH");
        AtomType typeHC = AtomType.element(Hydrogen.INSTANCE, "HC");
        AtomType typeHO = AtomType.element(Hydrogen.INSTANCE, "HO");
        species = new SpeciesBuilder(Space3D.getInstance())
                .addCount(typeCT135, 1)
                .addCount(typeCT157, 1)
                .addCount(typeOH, 1)
                .addCount(typeHC, 5)
                .addCount(typeHO, 1)
                .withConformation(conf)
                .build();
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularPeriodic(space, 100));

        SpeciesManager sm = getSpeciesManager();
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);

        double f = 0.01;
        Unit kJpmol = new UnitRatio(new PrefixedUnit(Prefix.KILO, Joule.UNIT), Mole.UNIT);
        P2Harmonic p2Bond = new P2Harmonic(kJpmol.toSim(f*224262.4), 1.529);
        List<int[]> pairs = new ArrayList<>();
        pairs.add(new int[]{0,1});
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        p2Bond = new P2Harmonic(kJpmol.toSim(f*267776), 1.41);
        pairs.clear();
        pairs.add(new int[]{1,2});
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        p2Bond = new P2Harmonic(kJpmol.toSim(f*284512), 1.09);
        pairs.clear();
        pairs.add(new int[]{0,3});
        pairs.add(new int[]{0,4});
        pairs.add(new int[]{0,5});
        pairs.add(new int[]{1,6});
        pairs.add(new int[]{1,7});
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        p2Bond = new P2Harmonic(kJpmol.toSim(f*462750.4), 0.945);
        pairs.clear();
        pairs.add(new int[]{2,8});
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        P3BondAngle p3Bond = new P3BondAngle(Math.PI*107.8/180.0, kJpmol.toSim(276.144));
        List<int[]> triplets = new ArrayList<>();
        triplets.add(new int[]{3,0,4});
        triplets.add(new int[]{3,0,5});
        triplets.add(new int[]{4,0,5});
        triplets.add(new int[]{6,1,7});
        pmBonding.setBondingPotentialTriplet(species, p3Bond, triplets);

        p3Bond = new P3BondAngle(Math.PI*110.7/180.0, kJpmol.toSim(313.8));
        triplets = new ArrayList<>();
        triplets.add(new int[]{1,0,3});
        triplets.add(new int[]{1,0,4});
        triplets.add(new int[]{1,0,5});
        triplets.add(new int[]{0,1,6});
        triplets.add(new int[]{0,1,7});
        pmBonding.setBondingPotentialTriplet(species, p3Bond, triplets);

        p3Bond = new P3BondAngle(Math.PI*109.5/180.0, kJpmol.toSim(418.4));
        triplets = new ArrayList<>();
        triplets.add(new int[]{0,1,2});
        pmBonding.setBondingPotentialTriplet(species, p3Bond, triplets);

        p3Bond = new P3BondAngle(Math.PI*109.5/180.0, kJpmol.toSim(292.88));
        triplets = new ArrayList<>();
        triplets.add(new int[]{6,1,2});
        triplets.add(new int[]{7,1,2});
        pmBonding.setBondingPotentialTriplet(species, p3Bond, triplets);

        p3Bond = new P3BondAngle(Math.PI*108.5/180.0, kJpmol.toSim(460.24));
        triplets = new ArrayList<>();
        triplets.add(new int[]{1,2,8});
        pmBonding.setBondingPotentialTriplet(species, p3Bond, triplets);

        P4BondTorsionPow p4Bond = new P4BondTorsionPow(new double[]{kJpmol.toSim(0.6276),
                                                                    kJpmol.toSim(1.8828),
                                                                    kJpmol.toSim(0),
                                                                    kJpmol.toSim(-2.5104)});
        List<int[]> quads = new ArrayList<>();
        quads.add(new int[]{3,0,1,6});
        quads.add(new int[]{4,0,1,6});
        quads.add(new int[]{5,0,1,6});
        quads.add(new int[]{3,0,1,7});
        quads.add(new int[]{4,0,1,7});
        quads.add(new int[]{5,0,1,7});
        pmBonding.setBondingPotentialQuad(species, p4Bond, quads);

        p4Bond = new P4BondTorsionPow(new double[]{kJpmol.toSim(0.97905),
                kJpmol.toSim(2.93716),
                kJpmol.toSim(0),
                kJpmol.toSim(-3.91622)});
        quads = new ArrayList<>();
        quads.add(new int[]{3,0,1,2});
        quads.add(new int[]{4,0,1,2});
        quads.add(new int[]{5,0,1,2});
        pmBonding.setBondingPotentialQuad(species, p4Bond, quads);

        p4Bond = new P4BondTorsionPow(new double[]{kJpmol.toSim(-0.4435),
                kJpmol.toSim(3.83255),
                kJpmol.toSim(0.72801),
                kJpmol.toSim(-4.11705)});
        quads = new ArrayList<>();
        quads.add(new int[]{0,1,2,8});
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

        TruncationFactory tf = new TruncationFactorySimple(rc);
        double[] epsilon = new double[5], sigma = new double[5];
        int idxCT135 = species.getTypeByName("CT135").getIndex();
        int idxCT157 = species.getTypeByName("CT157").getIndex();
        int idxOH = species.getTypeByName("OH").getIndex();
        int idxHO = species.getTypeByName("HO").getIndex();
        int idxHC = species.getTypeByName("HC").getIndex();
        sigma[idxCT135] = sigma[idxCT157] = 3.5;
        epsilon[idxCT135] = epsilon[idxCT157] = kJpmol.toSim(0.276144);
        sigma[idxOH] = 3.12;
        epsilon[idxOH] = kJpmol.toSim(0.71128);
        sigma[idxHC] = 2.5;
        epsilon[idxHC] = kJpmol.toSim(0.12552);
        List<P2LennardJones> allPLJ = new ArrayList<>();
        List<Double> allSigma = new ArrayList<>();
        for (int i=0; i<5; i++) {
            for (int j=i; j<5; j++) {
                double e = Math.sqrt(epsilon[i]*epsilon[j]);
                if (e == 0) continue;
                double s = (sigma[i]+sigma[j])/2;
                P2LennardJones p2LJ = new P2LennardJones(s, e);
                allPLJ.add(p2LJ);
                IPotential2 p2 = tf.make(p2LJ);
                allSigma.add(s);

                potentialMaster.setPairPotential(species.getAtomType(i), species.getAtomType(j), p2);
            }
        }

        potentialMaster.init();

        conf.minimize(box, pcAgg, species, 1e-10);

        box.setNMolecules(species, numMolecules);
        new BoxInflate(box, space, density).actionPerformed();

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
        potentialMaster.init();
        double u0 = potentialMaster.computeAll(false);
        double x = 1;
        while (u0 > 1e7*numMolecules) {
            x *= 0.99;
            for (int i=0; i<allPLJ.size(); i++) {
                allPLJ.get(i).setSigma(x*allSigma.get(i));
            }
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
                for (int i=0; i<allPLJ.size(); i++) {
                    allPLJ.get(i).setSigma(x*allSigma.get(i));
                }
                u0 = potentialMaster.computeAll(false);
            }
            integrator.reset();
        }
    }

    public static void main(String[] args) {

//        Unit kcalpmol = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
//        Unit kJpmol = new UnitRatio(new PrefixedUnit(Prefix.KILO, Joule.UNIT), Mole.UNIT);
//        double[] e = new double[]{0.276144, 0.12552, 0.71128};
//        System.out.println("epsilon");
//        for (int i=0; i<e.length; i++) {
//            System.out.println(e[i]+" "+kcalpmol.fromSim(kJpmol.toSim(e[i])));
//        }
//
//        System.out.println("bond angles");
//        double[] rd = new double[]{107.8, 110.7, 109.5, 108.5};
//        double[] r = new double[]{1.88146493365,1.93207948196,1.91113553093,1.89368223841};
//        for (int i=0; i<r.length; i++) {
//            System.out.println(r[i]+" "+Degree.UNIT.fromSim(r[i])+"   "+rd[i]+" "+Degree.UNIT.toSim(rd[i]));
//        }
//        System.exit(0);

        EthanolParams params = new EthanolParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.numSteps = 350000;
            params.graphics = true;
            params.numMolecules = 10;
            params.density = 0.01;
        }

        double mass = 2*Carbon.INSTANCE.getMass() + 6*Hydrogen.INSTANCE.getMass() + Oxygen.INSTANCE.getMass();
        Unit dUnit = new SimpleUnit(Null.DIMENSION, 1/(mass/Constants.AVOGADRO*1e24), "Density", "g/cm^3", false);

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

        final EthanolNPTMC sim = new EthanolNPTMC(Space3D.getInstance(), density, 5, numMolecules, temperature, rc, pressure);

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

            DiameterHashByElement dhbe = new DiameterHashByElement();
            dhbe.setDiameter("C", 3.5);
            dhbe.setDiameter("H", 2.5);
            dhbe.setDiameter("O", 3.12);
            simGraphic.getDisplayBox(sim.box).setDiameterHash(dhbe);
            ColorSchemeByType colorScheme = (ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme();
            colorScheme.setColor(sim.species.getTypeByName("CT135"), Color.GRAY);
            colorScheme.setColor(sim.species.getTypeByName("CT157"), Color.BLUE);
            colorScheme.setColor(sim.species.getTypeByName("OH"), Color.RED);
            colorScheme.setColor(sim.species.getTypeByName("HC"), Color.WHITE);
            colorScheme.setColor(sim.species.getTypeByName("HO"), Color.GREEN);

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

    public static class EthanolParams extends ParameterBase {
        public double temperatureK = 300;
        public int numMolecules = 500;
        public double density = 0.7893;
        public boolean graphics = false;
        public long numSteps = 200000;
        public double rc = 14;
        public double pressureKPa = 101;
    }
}
