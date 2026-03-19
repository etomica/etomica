package etomica.GasMOP;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDouble;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorListenerNHC;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveMoleculeRotate;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.function.FunctionMultiDimensionalDifferentiable;
import etomica.math.numerical.SteepestDescent;
import etomica.molecule.IMolecule;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.potential.UFF.*;
import etomica.potential.UFF.PDBReader;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeEwaldFourier;
import etomica.simulation.Simulation;
import etomica.simulation.prototypes.MCMoveWiggle;
import etomica.simulation.prototypes.MeterTorsionAngle;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.List;

public class GOMD  extends Simulation {
    public PotentialCompute pcAgg;
    public IntegratorMC integratorMC;
    public Box box;
    public MCMoveVolume mcMoveVolume;
    public IntegratorVelocityVerlet integrator;
    public MCMoveMolecule mcMoveMolecule;
    public ISpecies species;
    public PotentialMaster potentialMaster;
    double molecularWeight=0;
    public IntegratorListenerNHC nhc;
    public GOMD(Space space,int numMoleules, double temperature, String configFileName, Vector vecGrapheneone, Vector vecGraphenetwo, double x, double y, double z) {
        super(space);
        GrapheneReaderXYZPDB grapheneReaderXYZPDB = new GrapheneReaderXYZPDB();
        GOMDParams gomdParams = new GOMDParams();

        species = grapheneReaderXYZPDB.getSpecies(configFileName, new Vector3D(0, 0, 0), false, gomdParams.ifAmberData);
        System.out.println("Species");

        SpeciesBuilder speciesBuilder = new SpeciesBuilder(Space3D.getInstance());
        speciesBuilder.setDynamic(true);

     //   setRandom(new RandomMersenneTwister(3));
        addSpecies(species);

        box = this.makeBox();
        box.getBoundary().setBoxSize(new Vector3D(x,y,z)); //656050
        box.setNMolecules(species, 1);
        molecularWeight = species.getMass()*numMoleules;
        Vector vecZero = new Vector3D(0, 0, 10);
        Vector vecOne = new Vector3D(0, 0, 20);
        Vector vecTwo = new Vector3D(0, 0, -10);
        Vector vecThree = new Vector3D(0, 0, -20);
     /*   if (numMoleules > 1){
            System.out.println(vecZero);
            List<Vector> oldPositions = new ArrayList<>();
            IMolecule moleculeMOPZero = box.getMoleculeList().get(0);
            while (oldPositions.size() < moleculeMOPZero.getChildList().size()) {
                oldPositions.add(space.makeVector());
            }
            moleculeMOPZero.getChildList().forEach(atom -> {
                oldPositions.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(vecZero);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
            System.out.println(vecOne);
            IMolecule moleculeMOPOne = box.getMoleculeList().get(1);
            moleculeMOPOne.getChildList().forEach(atom -> {
                oldPositions.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(vecOne);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

            IMolecule moleculeMOPThree = box.getMoleculeList().get(2);
            moleculeMOPThree.getChildList().forEach(atom -> {
                oldPositions.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(vecTwo);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });

            System.out.println(vecGraphenetwo);
            IMolecule moleculeMOPFour = box.getMoleculeList().get(3);
            moleculeMOPFour.getChildList().forEach(atom -> {
                oldPositions.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(vecThree);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
        }*/

        double nbrRange = 12.5 * 1.05 + 1;
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);
        if(gomdParams.ifAmberData){
            grapheneReaderXYZPDB.makeBondingPotential(grapheneReaderXYZPDB, species, pmBonding, gomdParams.ifAmberData);
        } else {
            grapheneReaderXYZPDB.makeBondingPotential(grapheneReaderXYZPDB, species, pmBonding);
        }

        makeAtomPotentials(sm);
       // PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 5, pmBonding.getBondingInfo());
        potentialMaster = new PotentialMasterList(getSpeciesManager(), box, 2, nbrRange, pmBonding.getBondingInfo());
        SetPotential setPotential = new SetPotential();
        List<List<AtomType>> atomTypesGO =  setPotential.listFinal(species.getUniqueAtomTypes());
        //grapheneReaderXYZPDB.makeNBPotential(grapheneReaderXYZPDB, atomTypesGO, potentialMaster);

        PotentialComputeEwaldFourier ewaldFourier = new PotentialComputeEwaldFourier(getSpeciesManager(), box);
        PotentialComputeEwaldFourier.EwaldParams params = ewaldFourier.getOptimalParams(3, 0);

        List<AtomType> atomTypesMolecules = species.getAtomTypes();
        grapheneReaderXYZPDB.makeNBElectroPotential(grapheneReaderXYZPDB, atomTypesGO, potentialMaster, ewaldFourier, atomTypesMolecules, params);
        //System.exit(1);
        potentialMaster.doAllTruncationCorrection = true;
     //   pcAgg = new PotentialComputeAggregate(pmBonding, potentialMaster, ewaldFourier);
        pcAgg = new PotentialComputeAggregate(pmBonding, potentialMaster);
      //  pcAgg = new PotentialComputeAggregate(pmBonding);
        integrator = new IntegratorVelocityVerlet(pcAgg, random, 0.001, temperature, box);
        integrator.setIsothermal(true);
        integrator.setThermostatInterval(1000);
        integrator.setThermostatNoDrift(false);
        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
        nhc = new IntegratorListenerNHC(integrator, random, 3, 2);
        integrator.getEventManager().addListener(nhc);

       /* pcAgg = new PotentialComputeAggregate(pmBonding, potentialMaster);
        integratorMC = new IntegratorMC(pcAgg, random, temperature, box);
        getController().addActivity(new ActivityIntegrate(integrator));

        MCMoveAtom moveAtom = new MCMoveAtom(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(moveAtom);*/

        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
        potentialMaster.init();

      //  double u0 = potentialMasterCell.computeAll(false);
    /*    System.out.println("Before SD pmc : " + kcals.fromSim(u0));
        System.out.println("Before SD pmBonding: " + kcals.fromSim(u1));
        runSteepestDescentMinimization2(box, potentialMasterCell);
        u0 = potentialMasterCell.computeAll(false);
      //  u1 = pmBonding.computeAll(false);
        System.out.println("After Newton minimization: U pmc = " + kcals.fromSim(u0));
        System.out.println("After Newton minimization: U pmBonding = " + kcals.fromSim(u1));
        try{
            IMolecule moleculeGO = box.getMoleculeList().get(0);
            Map<Integer, Vector> positionMap = new HashMap<>();
            Map<Integer, String> atomMap = new HashMap<>();
            for (int i = 0; i < moleculeGO.getChildList().size(); i++){
                Vector vec = moleculeGO.getChildList().get(i).getPosition();
                AtomType atomName = moleculeGO.getChildList().get(i).getType();
                String atomTypeStringOne = String.valueOf(atomName);
                String stringAtomOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
                positionMap.put(i, vec);
                atomMap.put(i, stringAtomOne);
            }
            System.out.println(atomMap);
            System.out.println(positionMap);
            XYZWriter.writeXYZ(
                    "outMinimization2020.xyz",
                    atomMap,
                    positionMap,
                    "Generated from Etomica maps"
            );
        } catch (IOException e) {
            System.out.println("An error occurred while writing to the file: " + e.getMessage());
        }*/
        double uBonding = pmBonding.computeAll(false);
        double uAll = pcAgg.computeAll(false);
        System.out.println("utotal "+ uAll + " uBonding " + uBonding);
    }

    public static void main(String[] args) throws IOException {
        final String APP_NAME = "Methane Universal";
        GOMDParams params = new GOMDParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        Space space1 = Space3D.getInstance();
        Vector vector1 = new Vector3D(0, 0, 10);
        Vector vector2 = new Vector3D(0, 0, -10);
        Unit dUnit = new SimpleUnit(Null.DIMENSION, 1/(16.042/ Constants.AVOGADRO*1e24), "Density", "g/cm^3", false);

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

        GOMD sim = new GOMD(space1, numMolecules, temperature, configFilename, vector1, vector2, params.x, params.y, params.z );
        System.out.println("done writing ");
      //  System.exit(1);
        MeterPotentialEnergyFromIntegrator meterU = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        MeterTemperature meterT = new MeterTemperature(sim.box, 3);
        sim.potentialMaster.init();
        System.out.println("u0/N "+(meterU.getDataAsScalar()/numMolecules));
        Unit kjmol = new UnitRatio(new PrefixedUnit(Prefix.KILO,Joule.UNIT), Mole.UNIT);
        System.out.println("u0/N  "+ kjmol.fromSim(meterU.getDataAsScalar() / numMolecules) + " kJ/mol");
        // System.exit(1);

        MeterPressure meterP = new MeterPressure(sim.box, sim.pcAgg);
        meterP.setTemperature(temperature);
        meterP.doCallComputeAll(true);
        DataProcessorForked dpZ = new DataProcessorForked() {
            final DataDouble.DataInfoDouble dataInfo = new DataDouble.DataInfoDouble("Z", Null.DIMENSION);
            final DataDouble data = new DataDouble();

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
        System.out.println("Reached before dataforked");
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

        if (true) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Octane MD", 3);
            DiameterHashByType dhbt = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
            dhbt.setDiameter(sim.species.getAtomType(0), 3.75);
            dhbt.setDiameter(sim.species.getAtomType(1), 3.95);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            DataSourceCountTime timer = new DataSourceCountTime(sim.integrator);
            DisplayTextBox timerBox = new DisplayTextBox();
            timerBox.setLabel("Time");
            DataPumpListener pumpSteps = new DataPumpListener(timer, timerBox, 100);
            sim.integrator.getEventManager().addListener(pumpSteps);
            simGraphic.add(timerBox);

            Unit perN = new SimpleUnit(Null.DIMENSION, numMolecules, "1/N", "1/N", false);

            AccumulatorHistory historyU = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyU.setTimeDataSource(timer);
            AccumulatorHistory historyU2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyU2.setTimeDataSource(timer);
            AccumulatorAverageCollapsing avgEnergy = new AccumulatorAverageCollapsing();
            avgEnergy.setPushInterval(10);
            DataFork forkU = new DataFork(new IDataSink[]{historyU, historyU2, avgEnergy});
            DataPumpListener pumpU = new DataPumpListener(meterU, forkU, 10);
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

            if (sim.nhc != null) {
                IntegratorListenerNHC.DataSourceTotalEnergy meterTotalEnergy = new IntegratorListenerNHC.DataSourceTotalEnergy(sim.integrator, sim.nhc);
                AccumulatorHistory historyTotalEnergy = new AccumulatorHistory(new HistoryCollapsingDiscard());
                historyTotalEnergy.setTimeDataSource(timer);
                DataPumpListener pumpTotalEnergy = new DataPumpListener(meterTotalEnergy, historyTotalEnergy);
                sim.integrator.getEventManager().addListener(pumpTotalEnergy);
                DisplayPlotXChart plotTotalEnergy = new DisplayPlotXChart();
                plotTotalEnergy.setLabel("conserved energy");
                historyTotalEnergy.addDataSink(plotTotalEnergy.makeSink("total"));
                simGraphic.add(plotTotalEnergy);
            }

            AccumulatorHistory historyT = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyT.setTimeDataSource(timer);
            DataPumpListener pumpT = new DataPumpListener(meterT, historyT);
            sim.integrator.getEventManager().addListener(pumpT);
            DisplayPlotXChart plotT = new DisplayPlotXChart();
            plotT.setLabel("T");
            plotT.setUnit(Kelvin.UNIT);
            historyT.addDataSink(plotT.makeSink("T"));
            simGraphic.add(plotT);

            AccumulatorAverageCollapsing avgP = new AccumulatorAverageCollapsing();
            AccumulatorHistory historyP = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyP.setTimeDataSource(timer);
            AccumulatorHistory historyP2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyP2.setTimeDataSource(timer);
            forkP.addDataSink(avgP);
            forkP.addDataSink(historyP);
            forkP.addDataSink(historyP2);
            DataPumpListener pumpP = new DataPumpListener(meterP, forkP, 10);
            sim.integrator.getEventManager().addListener(pumpP);
            DisplayPlotXChart plotP = new DisplayPlotXChart();
            plotP.setLabel("P");
            historyP.addDataSink(plotP.makeSink("P"));
            plotP.setLegend(new DataTag[]{historyP.getTag()}, "samples");
            historyP2.addDataSink(plotP.makeSink("Pavg"));
            plotP.setLegend(new DataTag[]{historyP2.getTag()}, "avg");
            simGraphic.add(plotP);
            simGraphic.getController().getDataStreamPumps().add(pumpP);

            DisplayTextBoxesCAE displayP = new DisplayTextBoxesCAE();
            displayP.setAccumulator(avgP);
            simGraphic.add(displayP);

          /*  AccumulatorHistory historyZ = new AccumulatorHistory(new HistoryCollapsingDiscard());
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
            simGraphic.add(displayZ_);*/

            simGraphic.makeAndDisplayFrame();
            return;
        }
        File file = new File(params.fileOne);
        FileWriter writer = new FileWriter(file);

        System.out.println("Reached after for loop");

// Small equilibration
    /*    sim.getController().runActivityBlocking(
                new ActivityIntegrate(sim.integrator, numSteps/10)
        );*/

// Time source
        DataSourceCountTime timer = new DataSourceCountTime(sim.integrator);

// Energy history + averages
        AccumulatorHistory historyU = new AccumulatorHistory(new HistoryCollapsingDiscard());
        historyU.setTimeDataSource(timer);

        AccumulatorHistory historyU2 = new AccumulatorHistory(new HistoryCollapsingAverage());
        historyU2.setTimeDataSource(timer);

        AccumulatorAverageCollapsing avgEnergy = new AccumulatorAverageCollapsing();
        avgEnergy.setPushInterval(10);

        long samples = numSteps / (numMolecules* 8L);
        long bs = samples / 10;
        if (bs == 0) bs = 1;

        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(bs);
        forkP.addDataSink(accP);
        AccumulatorAverageFixed accZ = new AccumulatorAverageFixed(bs);
        dpZ.addDataSink(accZ);
        AccumulatorAverageFixed accZm1oR = new AccumulatorAverageFixed(bs);
        dpZm1oR.addDataSink(accZm1oR);
        DataPumpListener pumpP = new DataPumpListener(meterP, forkP, 8*numMolecules);
        sim.integrator.getEventManager().addListener(pumpP);

// Final block average
        AccumulatorAverageFixed accU =
                new AccumulatorAverageFixed((numSteps / 10) / 500);

// ⭐ ENERGY FILE WRITER
        EnergyHistoryWriter energyWriter =
                new EnergyHistoryWriter(writer, timer);

// Fork everything from meterU
        DataFork forkU = new DataFork(new IDataSink[]{
                historyU,
                historyU2,
                avgEnergy,
                accU,
                energyWriter
        });

// Pump every 10 steps
        DataPumpListener pumpU =
                new DataPumpListener(meterU, forkU, 10);

        sim.integrator.getEventManager().addListener(pumpU);

// ---- Run Production ----
        long t1 = System.nanoTime();

        sim.getController().runActivityBlocking(
                new ActivityIntegrate(sim.integrator, numSteps)
        );

        long t2 = System.nanoTime();

// Close writer
        energyWriter.close();
        writer.close();
        //sim.integratorMC.getMoveManager().addMCMove(sim.mcMoveMolecule);

        //sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps/5));





     //   AccumulatorAverageFixed accU = new AccumulatorAverageFixed((numSteps/10)/500);
     //  DataPumpListener pumpU = new DataPumpListener(meterU, accU, 200);
      //  sim.integrator.getEventManager().addListener(pumpU);



     /*   long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));
        long t2 = System.nanoTime();*/

        IData dataU = accU.getData();
        double avgU = dataU.getValue(AccumulatorAverage.AVERAGE.index) / numMolecules;
        double errU = dataU.getValue(AccumulatorAverage.ERROR.index) / numMolecules;
        double corU = dataU.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("U: " +" "+avgU+ "   err: "+" "+errU+"   cor: " +" "+corU);

        UnitRatio jouleMole = new UnitRatio(Joule.UNIT, Mole.UNIT);
        double valjouleMole = jouleMole.fromSim(dataU.getValue(AccumulatorAverage.AVERAGE.index/numMolecules));
        System.out.println(valjouleMole  / 1000 + " kJ/mol");
        IData dataP = accP.getData();
        UnitRatio den = new UnitRatio(Mole.UNIT, Liter.UNIT);
        System.out.println(den.fromSim(density) + " desnity");
        double avgP = dataP.getValue(AccumulatorAverage.AVERAGE.index);
        System.out.println(Bar.UNIT.fromSim(dataP.getValue(AccumulatorAverage.AVERAGE.index)) + " Pressure in Bar");
        double errP = dataP.getValue(AccumulatorAverage.ERROR.index);
        double corP = dataP.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("P: " +" "+Bar.UNIT.fromSim(avgP)+ "   err: " +" "+Bar.UNIT.fromSim(errP)+ "   cor: " +corP);

        IData dataZ = accZ.getData();
        double avgZ = dataZ.getValue(AccumulatorAverage.AVERAGE.index);
        double errZ = dataZ.getValue(AccumulatorAverage.ERROR.index);
        double corZ = dataZ.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("Z: "+" "+avgZ+"   err: "+" "+errZ+"   cor: "+" "+corZ);
        GOMDParams paramsN = new GOMDParams();
        IData dataZ_ = accZm1oR.getData();
        double avgZ_ = dataZ_.getValue(AccumulatorAverage.AVERAGE.index);
        double errZ_ = dataZ_.getValue(AccumulatorAverage.ERROR.index);
        double corZ_ = dataZ_.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);

        System.out.println("(Z-1)/rho: "+" "+avgZ_+"   err: "+" "+errZ_+"   cor: "+" "+corZ_);
        File file2 = new File(params.fileTwo);
        FileWriter writer2 = new FileWriter(file2);

            for (int m = 0; m < sim.box.getMoleculeList().size(); m++){
                IMolecule moleculeGO = sim.box.getMoleculeList().get(m);
                Map<Integer, Vector> positionMap = new HashMap<>();
                Map<Integer, String> atomMap = new HashMap<>();
                for (int i = 0; i < moleculeGO.getChildList().size(); i++){
                    Vector vec = moleculeGO.getChildList().get(i).getPosition();
                    AtomType atomName = moleculeGO.getChildList().get(i).getType();
                    String atomTypeStringOne = String.valueOf(atomName);
                    String stringAtomOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
                    positionMap.put(i, vec);
                    atomMap.put(i, stringAtomOne);
                }
                System.out.println(positionMap);
                System.out.println(atomMap);
                System.out.println("\n");
            }


          /*  System.out.println(atomMap);
            System.out.println(positionMap);
            XYZWriter.writeXYZ(
                   String.valueOf(file2),
                    atomMap,
                    positionMap,
                    "Generated from Etomica maps"
            );*/

        String absolutePath = file.getAbsolutePath();
        System.out.println("File path: " + absolutePath);
        System.out.println("time: "+" "+(t2-t1)/1e9);
    }

    public static class GOMDParams extends ParameterBase {
        public double temperatureK = 300;
        public int numMolecules = 1;
        //public int pressure = 10;
        public double density = 0.0000005;
        public boolean graphics = true;
        public long numSteps = 10000;
        public boolean ifAmberData = true;
        public double x =60;
        public double y =60;
        public double z = 50;
       // public String configFilename = "D:\\Sem-X\\GO\\graphitis\\GO_sheet2020";
        public String configFilename = "D:\\Sem-X\\GO\\amber\\go";
 //    public String configFilename = "GO_sheet";
     public String fileOne = "o1.txt";
     public String fileTwo = "o2.txt";
     //   public String outputFile = "D:\\Sem-X\\GO\\graphitis\\temp300K001Time.xyz";
        public int rc = 10;
        public double pressureKPa = 1402;
    }


    public static IPotential2[][] makeAtomPotentials(SpeciesManager sm) {
        // we could try to store the potentials more compactly, but it doesn't really matter
        ISpecies species = sm.getSpecies(sm.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        System.out.println(lastTypeIndex + 1+ " "+lastTypeIndex + 1 + " lastTypeIndex" + " "+species.getAtomType(species.getUniqueAtomTypeCount() - 1));
        return new IPotential2[lastTypeIndex + 1][lastTypeIndex + 1];
    }

    public void runSteepestDescentMinimization2(final Box box, final PotentialMasterCell pmc) {
        final IAtomList atoms = box.getLeafList();
        final int nAtoms = atoms.size();
        final int n = 3 * nAtoms;

        // 1) x0 from atom positions
        final double[] x0 = new double[n];
        for (int i = 0; i < nAtoms; i++) {
            Vector r = atoms.get(i).getPosition();
            x0[3*i]     = r.getX(0);
            x0[3*i + 1] = r.getX(1);
            x0[3*i + 2] = r.getX(2);
        }

        // 2) force buffer (Etomica-style). We'll reuse it to avoid allocations.
        final Vector[] forces = new Vector[nAtoms];
        for (int i = 0; i < nAtoms; i++) {
            forces[i] = box.getSpace().makeVector(); // or new Vector3D() depending on your space impl
        }

        // 3) Define f(x) and gradf(x)
        FunctionMultiDimensionalDifferentiable f = new FunctionMultiDimensionalDifferentiable() {

            @Override
            public int getDimension() {
                return n;
            }

            private void setPositionsFromX(double[] x) {
                for (int i = 0; i < nAtoms; i++) {
                    Vector r = atoms.get(i).getPosition();
                    r.setX(0, x[3*i]);
                    r.setX(1, x[3*i + 1]);
                    r.setX(2, x[3*i + 2]);
                }
            }

            @Override
            public double f(double[] x) {
                setPositionsFromX(x);
                // energy only
                return pmc.computeAll(false);
            }

            @Override
            public double[] gradf(double[] x) {
                final double h = 1e-6;
                double[] g = new double[x.length];
                double[] xx = x.clone();

                // IMPORTANT: make sure f(xx) updates atom positions from xx internally
                for (int i = 0; i < xx.length; i++) {
                    double old = xx[i];

                    xx[i] = old + h;
                    double fp = f(xx);

                    xx[i] = old - h;
                    double fm = f(xx);

                    xx[i] = old;
                    g[i] = (fp - fm) / (2*h);
                }
                return g;
            }

            // If your interface requires df(...) too, you can implement it
            // by delegating to gradf for order-1:
            @Override
            public double df(int[] d, double[] x) {
                int idx = -1, order = 0;
                for (int i = 0; i < d.length; i++) {
                    if (d[i] != 0) { order += d[i]; idx = i; }
                }
                if (order == 0) return f(x);
                if (order == 1 && idx >= 0) return gradf(x)[idx];
                throw new UnsupportedOperationException("Higher-order derivatives not supported");
            }
        };

        // 4) Run SD
        SteepestDescent sd = new SteepestDescent(f);

        double[] xStep = new double[n];
        Arrays.fill(xStep, 0.01);   // tune this; 1e-3 to 1e-2 is common for stiff bonded systems

        double tol = 1e-8;
        int maxIter = 500;

        double[] g0 = f.gradf(x0);
        double norm2 = 0;
        for (double gi : g0) norm2 += gi*gi;
        System.out.println("grad norm = " + Math.sqrt(norm2));
        double[] xmin = sd.minimize2(x0, xStep, tol, maxIter, true);
        double maxDx = 0;
        for (int i = 0; i < xmin.length; i++) {
            maxDx = Math.max(maxDx, Math.abs(xmin[i] - x0[i]));
        }
        System.out.println("max |dx| = " + maxDx);
        // 5) Copy back
        for (int i = 0; i < nAtoms; i++) {
            Vector r = atoms.get(i).getPosition();
            r.setX(0, xmin[3*i]);
            r.setX(1, xmin[3*i + 1]);
            r.setX(2, xmin[3*i + 2]);
        }
    }


}

