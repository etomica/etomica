package etomica.GasMOP;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataDouble;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveMoleculeRotate;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.potential.UFF.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.simulation.prototypes.MCMoveWiggle;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.VirialMultiUFF;


import java.util.*;
import java.util.List;

public class GasMOPSimulation extends Simulation {
    public PotentialCompute pcAgg;
    public IntegratorMC integratorMC;
    // public SpeciesGeneral species;
    public Box box;
    public MCMoveVolume mcMoveVolume;
    public MCMoveMolecule mcMoveMolecule;
    public ISpecies speciesMOP, speciesGas;
    double molecularWeight=0;
    protected Vector[] translationVectors;
    protected int[] constraintMap;
    public static int atom1, atom2, atom3, i=0;
    public static String atomName1, atomName2, atomName3;

    public GasMOPSimulation(Space space,String speciesMOPName,String speciesGasName, Vector centreMOP, Vector centreGas,  double density, int numMoleulesGas, int numMoleulesMOP, double temperature, String configFileName, int rc, double pressure){
        super(space);
        int i, j;

        //Call MoleculeOne
        System.out.println("Molecule One");

        speciesMOP = PDBReaderMOP.getSpecies(speciesMOPName);
        System.out.println(speciesMOP);
        printCOM(speciesMOP, centreMOP);
        List<AtomType> atomTypes1 = speciesMOP.getUniqueAtomTypes();
        List<List<AtomType>> pairsAtoms1 = VirialMultiUFF.getSpeciesPairs(speciesMOP);
        System.out.println(Arrays.deepToString(pairsAtoms1.toArray()));
        int pairAtomSize = pairsAtoms1.size();

        //Call Molecule Two
        System.out.println("Molecule Two");
        speciesGas = PDBReaderReplica.getSpecies(speciesGasName);
        printCOM(speciesGas, centreGas);
        List<List<AtomType>> pairsAtoms2 = VirialMultiUFF.getSpeciesPairs(speciesGas);
        int pairAtomSize2 = pairsAtoms2.size();
//System.exit(1);
        //Add both Molecules
        setRandom(new RandomMersenneTwister(1));
        addSpecies(speciesMOP);
        addSpecies(speciesGas);
        box = this.makeBox();
        box.setNMolecules(speciesMOP, numMoleulesMOP);
        box.setNMolecules(speciesGas, numMoleulesGas);
        new BoxInflate(box, space, density).actionPerformed();
        SpeciesManager sm1 = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).build();
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm1, box);
        UniversalSimulation.makeAtomPotentials(sm1);


        //Set Potential One
        ArrayList<ArrayList<Integer>> connectedAtoms1 = PDBReaderMOP.getConnectivityWithoutRunning();
        System.out.println(connectedAtoms1);
        ArrayList<ArrayList<Integer>> connectivityModified1 = PDBReaderMOP.getConnectivityModifiedWithoutRunning();
        System.out.println(connectivityModified1);
        Map<Integer,String> atomMap1 = PDBReaderMOP.getAtomMapWithoutRunning();
        System.out.println(atomMap1);
        HashMap<Integer, String> atomMapModified1 = PDBReaderMOP.getAtomMapModifiedWithoutRunning();
        System.out.println(atomMapModified1);
        ArrayList<Integer> bondList1 = PDBReaderMOP.getBondList(connectedAtoms1, atomMap1);
        System.out.println(bondList1);
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        Map<String, double[]> atomicPotMap1 = PDBReaderMOP.atomicPotMap();
        System.out.println(atomicPotMap1);
        ArrayList<Integer> bondsNum1 = PDBReaderMOP.getBonds();

        Map<Integer, String> atomIdentifierMapModified1 = PDBReaderMOP.atomIdentifierMapModified(connectivityModified1, atomMapModified1);
        List<int[]>dupletsSorted1= PDBReaderMOP.getDupletesSorted();
        List<int[]>tripletsSorted1= PDBReaderMOP.getAnglesSorted();
        List<int[]>quadrupletsSorted1= PDBReaderMOP.getTorsionSorted();

        Map<String[],List<int[]>> bondTypesMap1= PDBReaderMOP.idenBondTypes(dupletsSorted1, atomIdentifierMapModified1);
        Map<String[],List<int[]>> angleTypesMap1= PDBReaderMOP.idenAngleTypes(tripletsSorted1, atomIdentifierMapModified1);
        Map<String[],List<int[]>> torsionTypesMap1= PDBReaderMOP.idenTorsionTypes(quadrupletsSorted1, atomIdentifierMapModified1);
        ArrayList<ArrayList<Integer>> modifiedOutput1 = new ArrayList<>();

        for (ArrayList<Integer> innerList : connectivityModified1) {
            ArrayList<Integer> modifiedInnerList = new ArrayList<>(innerList.subList(1, innerList.size()));
            modifiedOutput1.add(modifiedInnerList);
        }
        System.out.println(modifiedOutput1 +" modified");
        IntArrayList[] dupletsIntArrayList1 = new IntArrayList[modifiedOutput1.size()];

        for (i = 0; i < modifiedOutput1.size(); i++) {
            ArrayList<Integer> innerList = modifiedOutput1.get(i);
            IntArrayList intArrayList = new IntArrayList(innerList.size());
            for (j = 0; j < innerList.size(); j++) {
                intArrayList.add(innerList.get(j));
               // System.out.println(intArrayList);
            }
            dupletsIntArrayList1[i] = intArrayList;
            //  System.out.println(dupletsIntArrayList[i]);
        }
        for (IntArrayList list : dupletsIntArrayList1) {
            for (i = 0; i < list.size(); i++) {
                int value = list.getInt(i);
                System.out.print(value + " ");
            }
        }
        PotentialMasterBonding.FullBondingInfo bondingInfo1 = new PotentialMasterBonding.FullBondingInfo(sm1);
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 5, pmBonding.getBondingInfo());
        doBondStrech(speciesMOP, bondTypesMap1, angleTypesMap1, torsionTypesMap1,bondsNum1,bondList1, quadrupletsSorted1, atomIdentifierMapModified1,atomicPotMap1, bondingInfo1);
        LJUFF[] p2LJ1 = new LJUFF[pairAtomSize];
        IPotential2[] p2lj1 = new IPotential2[pairAtomSize];
        double[] sigmaIJ1 = new double[pairAtomSize];
        doLJ(pairsAtoms1, p2LJ1, p2lj1, rc, potentialMasterCell, sigmaIJ1);


        //Set Potential Two
        ArrayList<ArrayList<Integer>> connectedAtoms2 =PDBReaderReplica.getConnectivityWithoutRunning();
        ArrayList<ArrayList<Integer>> connectivityModified2 = PDBReaderReplica.getConnectivityModifiedWithoutRunning();
        Map<Integer,String> atomMap2 = PDBReaderReplica.getAtomMapWithoutRunning();
        HashMap<Integer, String> atomMapModified2 = PDBReaderReplica.getAtomMapModifiedWithoutRunning();
        ArrayList<Integer> bondList2 = PDBReaderReplica.getBondList(connectedAtoms2, atomMap2);
        Map<String, double[]> atomicPotMap2 = PDBReaderReplica.atomicPotMap();
        Map<Integer, String> atomIdentifierMapModified2 = PDBReaderReplica.getatomIdentifierMapModified();

        List<int[]>dupletsSorted2= PDBReaderReplica.getDupletesSorted();
        List<int[]>tripletsSorted2=PDBReaderReplica.getAnglesSorted();
        List<int[]>quadrupletsSorted2=PDBReaderReplica.getTorsionSorted();
        ArrayList<Integer> bondsNum2 = PDBReaderReplica.getBonds();
        Map<String[],List<int[]>> bondTypesMap2= PDBReaderReplica.idenBondTypes(dupletsSorted2, atomIdentifierMapModified2);
        Map<String[],List<int[]>> angleTypesMap2= PDBReaderReplica.idenAngleTypes(tripletsSorted2, atomIdentifierMapModified2);
        Map<String[],List<int[]>> torsionTypesMap2= PDBReaderReplica.idenTorsionTypes(quadrupletsSorted2, atomIdentifierMapModified2);
        // System.out.println(connectivityModified2);

        ArrayList<ArrayList<Integer>> modifiedOutput2 = new ArrayList<>();
        for (ArrayList<Integer> innerList : connectivityModified2) {
            ArrayList<Integer> modifiedInnerList = new ArrayList<>(innerList.subList(1, innerList.size()));
            modifiedOutput2.add(modifiedInnerList);
        }
        System.out.println(modifiedOutput2 +" modified");
        IntArrayList[] dupletsIntArrayList2 = new IntArrayList[modifiedOutput2.size()];
        for (i = 0; i < modifiedOutput2.size(); i++) {
            ArrayList<Integer> innerList = modifiedOutput2.get(i);
            IntArrayList intArrayList = new IntArrayList(innerList.size());
            for (j = 0; j < innerList.size(); j++) {
                intArrayList.add(innerList.get(j));
               // System.out.println(intArrayList);
            }
            dupletsIntArrayList2[i] = intArrayList;
            //  System.out.println(dupletsIntArrayList[i]);
        }
        for (IntArrayList list : dupletsIntArrayList2) {
            for (i = 0; i < list.size(); i++) {
                int value = list.getInt(i);
                System.out.print(value + " ");
            }
        }
        doBondStrech(speciesGas,bondTypesMap2, angleTypesMap2, torsionTypesMap2, bondsNum2, bondList2,quadrupletsSorted2, atomIdentifierMapModified2, atomicPotMap2, bondingInfo1);
        LJUFF[] p2LJ2 = new LJUFF[pairAtomSize2];
        IPotential2[] p2lj2 = new IPotential2[pairAtomSize2];
        double[] sigmaIJ2 = new double[pairAtomSize2];
        doLJ(pairsAtoms2, p2LJ2, p2lj2, rc, potentialMasterCell, sigmaIJ2);


        //Set Interatomic Pot
        List<AtomType> list1 = speciesMOP.getUniqueAtomTypes();
        List<AtomType> list2 = speciesGas.getUniqueAtomTypes();
        List<AtomType> list3 = new ArrayList<>(list1);
        int list2Size = list2.size();
        boolean isEqual =false;
        for(i=0; i<list2Size; i++) {
            String name = list2.get(i).getName();
            for(j =0; j<list3.size(); j++){
                String nameSet = list3.get(j).getName();
                if(nameSet.equals(name)){
                    isEqual = true;
                    break;
                } else {
                    isEqual = false;
                }
            }
            if(!isEqual){
                list3.add(list2.get(i));
            }
        }
        System.out.println(list3);
        List<List<AtomType>> pairsAtomsTotal = new ArrayList<>();
        for(i=0; i<list3.size(); i++) {
            for (j = 0; j < list3.size(); j++) {
                if(i<=j){
                    List<AtomType> subPair = new ArrayList<>();
                    subPair.add(list3.get(i));
                    subPair.add(list3.get(j));
                    pairsAtomsTotal.add(subPair);
                }
            }
        }
        int pairAtomsTotalSize = pairsAtomsTotal.size();
        LJUFF[] p2LJTotal = new LJUFF[pairAtomsTotalSize];
        IPotential2[] p2ljTotal = new IPotential2[pairAtomsTotalSize];
        double[] sigmaIJTotal = new double[pairAtomsTotalSize];
        doLJ(pairsAtomsTotal, p2LJTotal, p2ljTotal, rc, potentialMasterCell, sigmaIJTotal);

        pcAgg = new PotentialComputeAggregate(pmBonding, potentialMasterCell);
        integratorMC = new IntegratorMC(pcAgg, random, temperature, box);
        getController().addActivity(new ActivityIntegrate(integratorMC));

        mcMoveMolecule = new MCMoveMolecule(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(mcMoveMolecule);

        MCMoveMoleculeRotate rotateMove = new MCMoveMoleculeRotate(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(rotateMove);

        MCMoveWiggle wiggleMove = new MCMoveWiggle(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(wiggleMove);

        MCMoveAtom moveAtom = new MCMoveAtom(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(moveAtom);


        if (configFileName != null) {
            ConfigurationFile config = new ConfigurationFile(configFileName);
            config.initializeCoordinates(box);
            BoxImposePbc.imposePBC(box);
        }
        else {
            ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
            configuration.initializeCoordinates(box);
            potentialMasterCell.init();
            double u0 = potentialMasterCell.computeAll(false);
            double x = 1;
           /* while (u0 > 1e4*numMoleulesGas){
                IMoleculeList moleculeList = box.getMoleculeList();
               // moleculeList.get(1).getChildList().forEach(atom -> {
                   // atom.getPosition().PE();
                });*/
                u0 = potentialMasterCell.computeAll(false);
                System.out.println(u0);

            integratorMC.reset();
            //System.out.println( u0 + " "+ x +" inMain "  + kcals.fromSim(u0));
            //System.out.println( u0 + " inMain afterwards "  + kcals.fromSim(u0));
        }
    }

    public static void main(String[] args) {
        final String APP_NAME = "Methane Universal";
        GasMOPSimulationParams params = new GasMOPSimulationParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.numSteps = 3500000;
            params.density = 0.00003;
            params.configFilename = null; // "octane";
            params.graphics = true;
        }


        Unit dUnit = new SimpleUnit(Null.DIMENSION, 1/(16.042/ Constants.AVOGADRO*1e24), "Density", "g/cm^3", false);

        double temperatureK = params.temperatureK;
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        System.out.println("Tsim "+temperature);
        int numMoleculesGas = params.numMoleculesGas;
        int numMoleculesMOP = params.NumMoleculesMOP;
        double density = dUnit.toSim(params.density);
        boolean graphics = params.graphics;
        long numSteps = params.numSteps;
        String configFilename = params.configFilename;
        int rc = params.rc;
        double pressureKPa = params.pressureKPa;
        Unit pUnit = new PrefixedUnit(Prefix.KILO, Pascal.UNIT);
        double pressure = pUnit.toSim(pressureKPa);
        String speciesMOPName = params.speciesMOPName;
        String speciesGasName = params.speciesGasName;
        Vector centreMOP = params.centreMOP;
        Vector centreGas = params.centreGas;

        System.out.println(numSteps+" steps");
        System.out.println("rc: "+rc);
        System.out.println("pressure "+ pressureKPa);
        System.out.println("initial density "+ density);
        System.out.println("initial density (g/cm^3) "+ dUnit.fromSim(density));
        final GasMOPSimulation sim = new GasMOPSimulation(Space3D.getInstance(),speciesMOPName, speciesGasName, centreMOP, centreGas, density, numMoleculesGas, numMoleculesMOP, temperature, configFilename, rc, pressure );
        MeterPotentialEnergyFromIntegrator meterU = new MeterPotentialEnergyFromIntegrator(sim.integratorMC);
        sim.integratorMC.getPotentialCompute().init();
        sim.integratorMC.reset();
        System.out.println("u0/N "+(meterU.getDataAsScalar()/numMoleculesGas));
        Unit kjmol = new UnitRatio(new PrefixedUnit(Prefix.KILO,Joule.UNIT), Mole.UNIT);
        System.out.println("u0/N  "+ kjmol.fromSim(meterU.getDataAsScalar() / numMoleculesGas) + " kJ/mol");
        System.exit(1);

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

       /*if (false) {
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

        }*/
        System.out.println("Reached after for loop");
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps/10));
        //sim.integratorMC.getMoveManager().addMCMove(sim.mcMoveMolecule);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps/5));

        long samples = numSteps / (numMoleculesGas* 8L);
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
        DataPumpListener pumpP = new DataPumpListener(meterP, forkP, 8*numMoleculesGas);
        sim.integratorMC.getEventManager().addListener(pumpP);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps));
        long t2 = System.nanoTime();

        IData dataU = accU.getData();
        double avgU = dataU.getValue(AccumulatorAverage.AVERAGE.index) / numMoleculesGas;
        double errU = dataU.getValue(AccumulatorAverage.ERROR.index) / numMoleculesGas;
        double corU = dataU.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("U: " +" "+avgU+ "   err: "+" "+errU+"   cor: " +" "+corU);

        UnitRatio jouleMole = new UnitRatio(Joule.UNIT, Mole.UNIT);
        double valjouleMole = jouleMole.fromSim(dataU.getValue(AccumulatorAverage.AVERAGE.index/numMoleculesGas));
        System.out.println(valjouleMole + " J/mol");
        IData dataP = accP.getData();
        UnitRatio den = new UnitRatio(Mole.UNIT, Liter.UNIT);
        System.out.println(den.fromSim(density) + " desnity");
        double avgP = dataP.getValue(AccumulatorAverage.AVERAGE.index);
        System.out.println(Bar.UNIT.fromSim(dataP.getValue(AccumulatorAverage.AVERAGE.index)) + " Pressure in Bar");
        double errP = dataP.getValue(AccumulatorAverage.ERROR.index);
        double corP = dataP.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("P: " +" "+avgP+ "   err: " +" "+errP+ "   cor: " +corP);

        IData dataZ = accZ.getData();
        double avgZ = dataZ.getValue(AccumulatorAverage.AVERAGE.index);
        double errZ = dataZ.getValue(AccumulatorAverage.ERROR.index);
        double corZ = dataZ.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("Z: "+" "+avgZ+"   err: "+" "+errZ+"   cor: "+" "+corZ);

        IData dataZ_ = accZm1oR.getData();
        double avgZ_ = dataZ_.getValue(AccumulatorAverage.AVERAGE.index);
        double errZ_ = dataZ_.getValue(AccumulatorAverage.ERROR.index);
        double corZ_ = dataZ_.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("(Z-1)/rho: "+" "+avgZ_+"   err: "+" "+errZ_+"   cor: "+" "+corZ_);

        System.out.println("time: "+" "+(t2-t1)/1e9);
    }

  /* public static void translateMolecule(ISpecies species, Box box){
        IRandom random = null;
        species.
       moleculeList.get(i).getChildList().forEach(atom -> {
           atom.getPosition().PE(translationVectors[tv]);
       });
    }*/
    public static void printCOM(ISpecies species, Vector center){
        Space space = Space3D.getInstance();
        //Vector center = new Vector3D(100,100,100);
        Vector dr = Vector.d(center.getD());
        double massSum = 0;
        IMolecule molecule = species.makeMolecule();
        IAtomList children = molecule.getChildList();
        int nAtoms = children.size();
        for (int i = 0; i < nAtoms; i++) {
            IAtom a = children.get(i);
            //System.out.println(a.getPosition() + " "+ i);
            double mass = a.getType().getMass();
            if (massSum == 0) {
                center.PEa1Tv1(mass, a.getPosition());
            } else {
                // sum = sum + mass*((sum/n)+pbc(r - sum/n))
                dr.E(a.getPosition());
                center.PEa1Tv1(mass, dr);
            }
            massSum += mass;
        }
        center.TE(1.0 / massSum);
         System.out.println(center + " 2 out");
    }

    public static void doBondStrech(ISpecies species1,Map<String[],List<int[]>> bondTypesMap1, Map<String[],List<int[]>> angleTypesMap1,Map<String[],List<int[]>> torsionTypesMap1,ArrayList<Integer> bondsNum1,ArrayList<Integer> bondList1, List<int[]>quadrupletsSorted1, Map<Integer, String> atomIdentifierMapModified1,Map<String, double[]> atomicPotMap1, PotentialMasterBonding.FullBondingInfo bondingInfo1){
        double Vi =0, Vj =0, V=0, Vtrue=0,  type;
        int p;
        int i =0;
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        for (Map.Entry<String[], List<int[]>> entry : bondTypesMap1.entrySet()) {
            String[] bondType = entry.getKey();
            List<int[]> bonds = entry.getValue();
            // System.out.println(Arrays.toString(bondType) + ": " + Arrays.deepToString(bonds.toArray()));
            for(int[]bondIndividual: bonds){
                bonds.add(bondIndividual);
                double[] bondParamsArray = new double[2];
                double[] bondConstant = new double[2];
                atom1 = bondIndividual[0];
                atom2 = bondIndividual[1];
                atomName1 = atomIdentifierMapModified1.get(atom1);
                atomName2 = atomIdentifierMapModified1.get(atom2);
                double[] atomOnePot = atomicPotMap1.get(atomName1);
                double[] atomTwoPot = atomicPotMap1.get(atomName2);
                double bondOrder = bondsNum1.get(i);
              /*  if(atomName1.equals("C_Ar") && atomName2.equals("C_Ar")){
                    bondOrder = 1.5;
                } else {
                    bondOrder = 1;
                }*/
                System.out.println(bondOrder +" bondorder");
                //  System.out.println(Arrays.toString(dupletsSorted.get(i)) + " " + bondOrder+ " " + atomName1 + " " + atomName2+" "+ Arrays.toString(atomOnePot) +" " + Arrays.toString(atomTwoPot));
                bondParamsArray= UFF.bondUFF (atomOnePot[0],  atomTwoPot[0], atomOnePot[5],  atomTwoPot[5], atomOnePot[6], atomTwoPot[6], bondOrder);
                //System.out.println(Arrays.toString(bondParamsArray) + " ArrayToString");
                bondConstant = UFF.BondConstantArray(bondParamsArray[0], bondParamsArray[1]);
                //System.out.println(Arrays.toString(bondConstant) + " arrayConstant");
                P2HarmonicUFF p2Bond = new P2HarmonicUFF(bondParamsArray[0],  bondParamsArray[1]);
                //bondParams.add(bondConstant);
                bondingInfo1.setBondingPotentialPair(species1, p2Bond, bonds);
                i++;
                break;
            }
        }

        for (Map.Entry<String[], List<int[]>> entry : angleTypesMap1.entrySet()) {
            String[] angleType = entry.getKey();
            List<int[]> angle = entry.getValue();
            // System.out.println(Arrays.toString(angleType) + ": " + Arrays.deepToString(angle.toArray()));
            for(int[]angleIndividual: angle){
                double[] angleParamsArray = new double[4];
                atom1 = angleIndividual[0];
                atom2 = angleIndividual[1];
                atom3 = angleIndividual[2];
                atomName1 = atomIdentifierMapModified1.get(atom1);
                atomName2 = atomIdentifierMapModified1.get(atom2);
                atomName3 = atomIdentifierMapModified1.get(atom3);
                int bondListValueOne = bondList1.get(atom1);
                int bondListValueTwo = bondList1.get(atom2);
                int bondListValueThree = bondList1.get(atom3);
                double[] atomOnePot = atomicPotMap1.get(atomName1);
                double[] atomTwoPot = atomicPotMap1.get(atomName2);
                double[] atomThreePot = atomicPotMap1.get(atomName3);
                int num =0;
                int caseNum =1;
                angleParamsArray= UFF.angleUFF (atomOnePot[0], atomTwoPot[0], atomThreePot[0], atomOnePot[5], atomTwoPot[5], atomThreePot[5], atomOnePot[6], atomTwoPot[6],atomThreePot[6], atomTwoPot[1], bondListValueOne, bondListValueTwo, bondListValueThree,0);
                // System.out.println(Arrays.toString(angleParamsArray) + " arrayAngle");
                P3BondAngleUFF p3Angle = new P3BondAngleUFF(angleParamsArray[0],  angleParamsArray[1], angleParamsArray[2], angleParamsArray[3],atomTwoPot[1], 0, caseNum);
                bondingInfo1.setBondingPotentialTriplet(species1, p3Angle, angle);
                break;
            }
        }

        P4BondTorsionUFF[] p4BondTorsionArray2 = new P4BondTorsionUFF[quadrupletsSorted1.size()];
        ArrayList<double[]> p4ValueArray2 = new ArrayList<>();

        for (Map.Entry<String[], List<int[]>> entry : torsionTypesMap1.entrySet()) {
            i = 0;
            String[] torsionType = entry.getKey();
            List<int[]> torsion = entry.getValue();
            // System.out.println(Arrays.toString(torsionType) + ": " + Arrays.deepToString(torsion.toArray()));
            for(int[]torsionIndividual: torsion){
                type = 0;
                p = 0;
                double[] torsionParamsArray = new double[4];
                atom2 = torsionIndividual[1];
                atom3 = torsionIndividual[2];
                atomName2 = atomIdentifierMapModified1.get(atom2);
                atomName3 = atomIdentifierMapModified1.get(atom3);
                Vi = UFF.switchCaseTorsion(atomName2);
                int bondListValueOne = bondList1.get(torsionIndividual[1]);
                int bondListValueTwo = bondList1.get(torsionIndividual[2]);
                p = p + bondListValueOne + bondListValueTwo;
                Vj = UFF.switchCaseTorsion(atomName3);
                V = Math.sqrt(Vi*Vj);
                Vtrue = kcals.toSim(V);
                double bondOrder = 1 ;
                torsionParamsArray = UFF.torsionUFF(Vtrue, p, bondOrder);
                p4BondTorsionArray2[i] = new P4BondTorsionUFF(torsionParamsArray[0], (int) torsionParamsArray[1], torsionParamsArray[2]);
                double[] array = {torsionParamsArray[0], torsionParamsArray[1], torsionParamsArray[2]};
                p4ValueArray2.add(array);
                //System.out.println(Arrays.toString(array));
                bondingInfo1.setBondingPotentialQuad(species1, p4BondTorsionArray2[i], torsion);
                i++;
            }
        }
    }
    public static void doLJ(List<List<AtomType>> pairsAtoms, LJUFF[] p2LJ, IPotential2[] p2lj, int rc, PotentialMasterCell potentialMasterCell, double[] sigmaIJ){
        int i = 0;
        UFF uff = new UFF();
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        System.out.println(pairsAtoms + " Pairs");
        for(List<AtomType>individualPair: pairsAtoms){
            AtomType atomNameOne = individualPair.get(0);
            AtomType atomNameTwo = individualPair.get(1);
            String atomTypeStringOne = String.valueOf(atomNameOne);
            String atomTypeStringTwo = String.valueOf(atomNameTwo);
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);
            double[] iKey = PDBReaderReplica.atomicPot(atomTypeOne);
            double[] jKey = PDBReaderReplica.atomicPot(atomTypeTwo);
            double epsilonIKey = kcals.toSim(iKey[3]);
            double epsilonJKey = kcals.toSim(jKey[3]);
            double sigmaIKey = iKey[2];
            double sigmaJKey = jKey[2];
            sigmaIJ[i] = (sigmaIKey + sigmaJKey) / 2;
            TruncationFactory tf = new TruncationFactoryForceShift(rc);
            p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);
            p2lj[i] = tf.make(p2LJ[i]);
            // potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo,p2lj[i]);
            potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2lj[i], new double[]{1, 0, 0, 1});
            i++;
        }
    }


    public static class GasMOPSimulationParams extends ParameterBase {
        public String speciesMOPName = "F://12";
        public String speciesGasName = "F://Avagadro//molecule//ch4";
        public Vector centreGas = new Vector3D(2000,2000,2000);
        public Vector centreMOP = new Vector3D(-2000,-2000,-2000);
        public double temperatureK = 300;
        public int numMoleculesGas = 1;
        public int NumMoleculesMOP = 1;
        //public int pressure = 10;
        public double density = 0.0000005;
        public boolean graphics = false;
        public long numSteps = 1;
        public String configFilename = null;
        public int rc = 10;
        public double pressureKPa = 1402;

    }
}
