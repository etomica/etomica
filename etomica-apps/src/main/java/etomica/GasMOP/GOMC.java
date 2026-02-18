package etomica.GasMOP;

import etomica.action.IAction;
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
import etomica.integrator.IntegratorListenerNHC;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.IMolecule;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.IPotential2;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
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
import etomica.virial.mcmove.MCMoveClusterAtomMulti;

import java.awt.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GOMC  extends Simulation {
    public PotentialCompute pcAgg;
    public IntegratorMC integratorMC;
    public Box box;
    public MCMoveVolume mcMoveVolume;
    public MCMoveMolecule mcMoveMolecule;
    public ISpecies species;
    double molecularWeight=0;
    public IntegratorListenerNHC nhc;
    public GOMC(Space space, int numMoleules, double temperature, String configFileName, Vector vecGrapheneone, Vector vecGraphenetwo) {
        super(space);
        GrapheneReaderXYZPDB grapheneReaderXYZPDB = new GrapheneReaderXYZPDB();
        species = grapheneReaderXYZPDB.getSpecies(configFileName, new Vector3D(0, 0, 0), false);

        System.out.println("Species");

        SpeciesBuilder speciesBuilder = new SpeciesBuilder(Space3D.getInstance());
        speciesBuilder.setDynamic(true);

        setRandom(new RandomMersenneTwister(1));
        addSpecies(species);

        box = this.makeBox();
        box.getBoundary().setBoxSize(new Vector3D(31,31,30));
        box.setNMolecules(species, numMoleules);
        molecularWeight = species.getMass()*numMoleules;
        if (numMoleules > 1){
            System.out.println(vecGrapheneone);
            List<Vector> oldPositions = new ArrayList<>();
            IMolecule moleculeMOPZero = box.getMoleculeList().get(0);
            while (oldPositions.size() < moleculeMOPZero.getChildList().size()) {
                oldPositions.add(space.makeVector());
            }
            moleculeMOPZero.getChildList().forEach(atom -> {
                oldPositions.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(vecGrapheneone);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
            System.out.println(vecGraphenetwo);
            IMolecule moleculeMOPOne = box.getMoleculeList().get(1);
            moleculeMOPOne.getChildList().forEach(atom -> {
                oldPositions.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(vecGraphenetwo);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
        }


        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);
        grapheneReaderXYZPDB.makeBondingPotential(grapheneReaderXYZPDB, species, pmBonding);
        makeAtomPotentials(sm);
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 5, pmBonding.getBondingInfo());
        SetPotential setPotential = new SetPotential();
        List<List<AtomType>> atomTypesGO =  setPotential.listFinal(species.getUniqueAtomTypes());
        grapheneReaderXYZPDB.makeNBPotential(grapheneReaderXYZPDB, atomTypesGO, potentialMasterCell);

        potentialMasterCell.doAllTruncationCorrection = false;


        pcAgg = new PotentialComputeAggregate(pmBonding, potentialMasterCell);
        integratorMC = new IntegratorMC(pcAgg, random, temperature, box);
        getController().addActivity(new ActivityIntegrate(integratorMC));

        MCMoveAtom moveAtom = new MCMoveAtom(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(moveAtom);
        moveAtom.setStepSize(5);

        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT),Mole.UNIT);
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
        potentialMasterCell.init();


        double u0 = potentialMasterCell.computeAll(false);
        System.out.println( u0 +" initial Value "  + kcals.fromSim(u0));

    }

    public static void main(String[] args) throws IOException {
        final String APP_NAME = "Methane Universal";
        GOMCParams params = new GOMCParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        Space space1 = Space3D.getInstance();
        Vector vector1 = new Vector3D(0, 0, 5);
        Vector vector2 = new Vector3D(0, 0, -5);
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

        GOMC sim = new GOMC(space1, 1, 500, configFilename, vector1, vector2  );
        System.out.println("done writing ");
        //  System.exit(1);
        MeterPotentialEnergyFromIntegrator meterU = new MeterPotentialEnergyFromIntegrator(sim.integratorMC);
        sim.integratorMC.getPotentialCompute().init();
        sim.integratorMC.reset();
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
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorMC, 5000000));
            sim.getController().addActionSequential(new IAction() {
                @Override
                public void actionPerformed() {
                    sim.integratorMC.getMoveManager().addMCMove(sim.mcMoveVolume);
                }
            });
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorMC));
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Universal MC", 3);
            System.out.println("Reached after simulation graphic");
            DiameterHashByType dhbt = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
            ((DiameterHashByType)((DisplayBox)simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getAtomType(0),1);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getAtomType(0), Color.darkGray);
            /*((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getTypeByName("C_1"), Color.lightGray);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getTypeByName("O_2"), Color.red);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getTypeByName("O_1"), ColorExtra.indianRed);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getTypeByName("H_2"), ColorExtra.cornflowerblue);*/

            DataSourceCountSteps timeSource = new DataSourceCountSteps(sim.integratorMC);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            List<DataPump> dataPumps = simGraphic.getController().getDataStreamPumps();
            // System.out.println("Reached dataPump");
            DataSourceCountSteps timer = new DataSourceCountSteps(sim.integratorMC);
            DisplayTextBox timerBox = new DisplayTextBox();
            timerBox.setLabel("Steps");
            DataPumpListener pumpSteps = new DataPumpListener(timer, timerBox, numMolecules);
            sim.integratorMC.getEventManager().addListener(pumpSteps);
            simGraphic.add(timerBox);

            MeterDensity meterDensity = new MeterDensity(sim.box());
            AccumulatorHistory accDensity = new AccumulatorHistory(new HistoryCollapsingAverage());
            accDensity.setTimeDataSource(timer);
            DataPumpListener pumpDensity = new DataPumpListener(meterDensity, accDensity, 10);
            sim.integratorMC.getEventManager().addListener(pumpDensity);
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
            sim.integratorMC.getEventManager().addListener(pumpU);
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
            sim.integratorMC.getEventManager().addListener(pumpP);
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
            simGraphic.makeAndDisplayFrame(APP_NAME);


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
            sim.integratorMC.getEventManager().addListener(pumpTorsion);
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
        File file = new File("output.txt");
        FileWriter writer = new FileWriter(file);
        System.out.println("Reached after for loop");
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps/10));
        //sim.integratorMC.getMoveManager().addMCMove(sim.mcMoveMolecule);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps/5));

        long samples = numSteps / (numMolecules* 8L);
        long bs = samples / 10;
        if (bs == 0) bs = 1;

        AccumulatorAverageFixed accU = new AccumulatorAverageFixed((numSteps/10)/500);
        DataPumpListener pumpU = new DataPumpListener(meterU, accU, 200);
        sim.integratorMC.getEventManager().addListener(pumpU);

        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(bs);
        forkP.addDataSink(accP);
        AccumulatorAverageFixed accZ = new AccumulatorAverageFixed(bs);
        dpZ.addDataSink(accZ);
        AccumulatorAverageFixed accZm1oR = new AccumulatorAverageFixed(bs);
        dpZm1oR.addDataSink(accZm1oR);
        DataPumpListener pumpP = new DataPumpListener(meterP, forkP, 8*numMolecules);
        sim.integratorMC.getEventManager().addListener(pumpP);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps));
        long t2 = System.nanoTime();

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
        GOMD.GOMDParams paramsN = new GOMD.GOMDParams();
        IData dataZ_ = accZm1oR.getData();
        double avgZ_ = dataZ_.getValue(AccumulatorAverage.AVERAGE.index);
        double errZ_ = dataZ_.getValue(AccumulatorAverage.ERROR.index);
        double corZ_ = dataZ_.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("(Z-1)/rho: "+" "+avgZ_+"   err: "+" "+errZ_+"   cor: "+" "+corZ_);
        try{
            IMolecule moleculeGO = sim.box.getMoleculeList().get(0);
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
                    paramsN.outputFile,
                    atomMap,
                    positionMap,
                    "Generated from Etomica maps"
            );
        } catch (IOException e) {
            System.out.println("An error occurred while writing to the file: " + e.getMessage());
        }
        writer.close();
        String absolutePath = file.getAbsolutePath();
        System.out.println("File path: " + absolutePath);
        System.out.println("time: "+" "+(t2-t1)/1e9);
    }

    public static class GOMCParams extends ParameterBase {
        public double temperatureK = 500;
        public int numMolecules = 1;
        //public int pressure = 10;
        public double density = 0.0000005;
        public boolean graphics = false;
        public long numSteps = 10000000;
        public String configFilename = "D:\\Sem-X\\GO\\graphitis\\GOsheet30";
        public String outputFile = "D:\\Sem-X\\GO\\graphitis\\temp30L5001M.xyz";
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
}
