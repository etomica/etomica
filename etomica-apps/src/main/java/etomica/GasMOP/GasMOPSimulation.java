package etomica.GasMOP;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.box.Box;
import etomica.box.RandomPositionSourceRectangular;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterWidomInsertion;
import etomica.data.types.DataDouble;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.*;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.function.IFunction;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSourceRandomMolecule;
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
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.VirialMultiUFF;


import java.awt.*;
import java.util.*;
import java.util.List;

public class GasMOPSimulation extends Simulation{
    public PotentialCompute pcAgg;
    public IntegratorMC integratorMC;
    public Box box;
    public MCMoveVolume mcMoveVolume;
    public MCMoveInsertDelete mcMoveID;
    public MCMoveMolecule mcMoveMolecule;
    public MCMoveMoleculeRotate mcMoveMoleculeRotate;
    public MCMoveWiggle wiggleMove;
   // public MCMoveAtom mcMoveAtom;
    public ISpecies speciesMOP, speciesGas;
    double molecularWeight=0;
    protected Vector[] translationVectors;
    protected int[] constraintMap;
    public static int atom1, atom2, atom3, i=0;
    public static String atomName1, atomName2, atomName3;


    public GasMOPSimulation(Space space,String speciesMOPName, Vector centreMOP,String speciesGasName,  double density, int numMoleulesGas, int numMoleulesMOP, double temperature, String configFileName, double rc, double pressure, boolean isDynamic, Vector boxSize){
        super(space);
        int i, j;
        System.out.println(rc);
       // setRandom(new RandomMersenneTwister(1));
        //Call MoleculeOne
  //      System.out.println("Molecule One");
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
        SetPotential doPotential = new SetPotential();
        speciesMOP = pdbReaderMOP.getSpecies(speciesMOPName, isDynamic, centreMOP, false, false);
      //  System.out.println(speciesMOP);
       // printCOM(speciesMOP, centreMOP);
        List<AtomType> atomTypes1 = speciesMOP.getUniqueAtomTypes();
        List<List<AtomType>> pairsAtoms1 = VirialMultiUFF.getSpeciesPairs(speciesMOP);
       // System.out.println(Arrays.deepToString(pairsAtoms1.toArray()));
        int pairAtomSize = pairsAtoms1.size();

        //Call Molecule Two
      //  System.out.println("Molecule Two");
        speciesGas = pdbReaderReplica.getSpecies(speciesGasName, true,new Vector3D(0,0,0), false);
     //   printCOM(speciesGas, centreGas);
        List<List<AtomType>> pairsAtoms2 = VirialMultiUFF.getSpeciesPairs(speciesGas);
        int pairAtomSize2 = pairsAtoms2.size();
        //Add both Molecules

        addSpecies(speciesMOP);
        addSpecies(speciesGas);
        box = this.makeBox();
        box.setNMolecules(speciesMOP, numMoleulesMOP);
        new BoxInflate(box, space, density).actionPerformed();

        SpeciesManager sm1 = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).build();
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm1, box);
        doPotential.makeAtomPotentials(sm1);


        //Set Potential One
        ArrayList<ArrayList<Integer>> connectedAtoms1 = pdbReaderMOP.getConnectivity();
        //System.out.println(connectedAtoms1);
        ArrayList<ArrayList<Integer>> connectivityModified1 = pdbReaderMOP.getConnectivityModifiedWithoutRunning();
       // System.out.println(connectivityModified1);
        Map<Integer,String> atomMap1 = pdbReaderMOP.getAtomMapWithoutRunning();
       // System.out.println(atomMap1);
        Map<Integer, String> atomMapModified1 = pdbReaderMOP.getAtomMapModifiedWithoutRunning();
       // System.out.println(atomMapModified1);
        ArrayList<Integer> bondList1 = pdbReaderMOP.getBondList(connectedAtoms1, atomMap1);
       // System.out.println(bondList1);
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        Map<String, double[]> atomicPotMap1 = pdbReaderMOP.atomicPotMap();
       // System.out.println(atomicPotMap1);
        ArrayList<Integer> bondsNum1 = pdbReaderMOP.getBonds();

        Map<Integer, String> atomIdentifierMapModified1 = pdbReaderMOP.getModifiedAtomIdentifierMap();
        List<int[]>dupletsSorted1= pdbReaderMOP.getDupletesSorted();
        List<int[]>tripletsSorted1= pdbReaderMOP.getAnglesSorted();
        List<int[]>quadrupletsSorted1= pdbReaderMOP.getTorsionSorted();

        Map<String[],List<int[]>> bondTypesMap1= pdbReaderMOP.idenBondTypes(dupletsSorted1, atomIdentifierMapModified1);
        Map<String[],List<int[]>> angleTypesMap1= pdbReaderMOP.idenAngleTypes(tripletsSorted1, atomIdentifierMapModified1);
        Map<String[],List<int[]>> torsionTypesMap1= pdbReaderMOP.idenTorsionTypes(quadrupletsSorted1, atomIdentifierMapModified1);
        ArrayList<ArrayList<Integer>> modifiedOutput1 = new ArrayList<>();

        for (ArrayList<Integer> innerList : connectivityModified1) {
            ArrayList<Integer> modifiedInnerList = new ArrayList<>(innerList.subList(1, innerList.size()));
            modifiedOutput1.add(modifiedInnerList);
        }
     //   System.out.println(modifiedOutput1 +" modified");
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
       /* for (IntArrayList list : dupletsIntArrayList1) {
            for (i = 0; i < list.size(); i++) {
                int value = list.getInt(i);
               // System.out.print(value + " ");
            }
        }*/
        SetPotential setPotential = new SetPotential();
        PotentialMasterBonding.FullBondingInfo bondingInfo1 = new PotentialMasterBonding.FullBondingInfo(sm1);
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 5, pmBonding.getBondingInfo());
        setPotential.setupBondStrech(speciesMOP, bondTypesMap1, angleTypesMap1, torsionTypesMap1,bondsNum1,bondList1, quadrupletsSorted1, atomIdentifierMapModified1,atomicPotMap1, bondingInfo1, pmBonding);
        LJUFF[] p2LJ1 = new LJUFF[pairAtomSize];
        IPotential2[] p2lj1 = new IPotential2[pairAtomSize];
        double[] sigmaIJ1 = new double[pairAtomSize];
        //doPotential.doLJ(pairsAtoms1, p2LJ1, p2lj1, rc, potentialMasterCell, sigmaIJ1);


        //Set Potential Two
        ArrayList<ArrayList<Integer>> connectedAtoms2 =pdbReaderReplica.getConnectivityWithoutRunning();
        ArrayList<ArrayList<Integer>> connectivityModified2 = pdbReaderReplica.getConnectivityModifiedWithoutRunning();
        Map<Integer,String> atomMap2 = pdbReaderReplica.getAtomMapWithoutRunning();
        HashMap<Integer, String> atomMapModified2 = pdbReaderReplica.getAtomMapModifiedWithoutRunning();
        ArrayList<Integer> bondList2 = pdbReaderReplica.getBondList(connectedAtoms2, atomMap2);
        Map<String, double[]> atomicPotMap2 = pdbReaderReplica.atomicPotMap();
        Map<Integer, String> atomIdentifierMapModified2 = pdbReaderReplica.getatomIdentifierMapModified();

        List<int[]>dupletsSorted2= pdbReaderReplica.getDupletesSorted();
        List<int[]>tripletsSorted2=pdbReaderReplica.getAnglesSorted();
        List<int[]>quadrupletsSorted2=pdbReaderReplica.getTorsionSorted();
        ArrayList<Integer> bondsNum2 = pdbReaderReplica.getBonds();
        Map<String[],List<int[]>> bondTypesMap2= pdbReaderReplica.idenBondTypes(dupletsSorted2, atomIdentifierMapModified2);
        Map<String[],List<int[]>> angleTypesMap2= pdbReaderReplica.idenAngleTypes(tripletsSorted2, atomIdentifierMapModified2);
        Map<String[],List<int[]>> torsionTypesMap2= pdbReaderReplica.idenTorsionTypes(quadrupletsSorted2, atomIdentifierMapModified2);
        // System.out.println(connectivityModified2);

        ArrayList<ArrayList<Integer>> modifiedOutput2 = new ArrayList<>();
        for (ArrayList<Integer> innerList : connectivityModified2) {
            ArrayList<Integer> modifiedInnerList = new ArrayList<>(innerList.subList(1, innerList.size()));
            modifiedOutput2.add(modifiedInnerList);
        }
     //   System.out.println(modifiedOutput2 +" modified");
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
      /*  for (IntArrayList list : dupletsIntArrayList2) {
            for (i = 0; i < list.size(); i++) {
                int value = list.getInt(i);
                System.out.print(value + " ");
            }
        }*/
        setPotential.setupBondStrech(speciesGas,bondTypesMap2, angleTypesMap2, torsionTypesMap2, bondsNum2, bondList2,quadrupletsSorted2, atomIdentifierMapModified2, atomicPotMap2, bondingInfo1, pmBonding);
        LJUFF[] p2LJ2 = new LJUFF[pairAtomSize2];
        IPotential2[] p2lj2 = new IPotential2[pairAtomSize2];
        double[] sigmaIJ2 = new double[pairAtomSize2];
        //doPotential.doLJ(pairsAtoms2, p2LJ2, p2lj2, rc, potentialMasterCell, sigmaIJ2);


        //Set Interatomic Pot
        List<AtomType> list1 = speciesMOP.getUniqueAtomTypes();
        List<AtomType> list2 = speciesGas.getUniqueAtomTypes();
        List<AtomType> list3 = new ArrayList<>(list1);
        int list2Size = list2.size();
        boolean isEqual =false;
        List<List<AtomType>> pairsAtomsTotal = new ArrayList<>();
        for( i=0; i<list1.size(); i++){
            for ( j=0; j<list2.size(); j++){
                List<AtomType> subPair = new ArrayList<>();
                subPair.add(list1.get(i));
                subPair.add(list2.get(j));
                pairsAtomsTotal.add(subPair);
            }
        }
      /*  for(i=0; i<list2Size; i++) {
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
        }*/
        //System.out.println(pairsAtomsTotal);
        int pairAtomsTotalSize = pairsAtomsTotal.size();
        LJUFF[] p2LJTotal = new LJUFF[pairAtomsTotalSize];
        IPotential2[] p2ljTotal = new IPotential2[pairAtomsTotalSize];
        double[] sigmaIJTotal = new double[pairAtomsTotalSize];
        doPotential.doLJ(pairsAtomsTotal, p2LJTotal, p2ljTotal, rc, potentialMasterCell, sigmaIJTotal);
        box.getBoundary().setBoxSize(boxSize);
        pcAgg = new PotentialComputeAggregate(pmBonding, potentialMasterCell);
        integratorMC = new IntegratorMC(pcAgg, random, temperature, box);
        getController().addActivity(new ActivityIntegrate(integratorMC));

        MCMoveAtom mcMoveAtom = new MCMoveAtom(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(mcMoveAtom);
        //((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesMOP);
        /*mcMoveMolecule = new MCMoveMolecule(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(mcMoveMolecule);
        ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesGas);

        mcMoveMoleculeRotate = new MCMoveMoleculeRotate(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(mcMoveMoleculeRotate);
        ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesGas);*/

        potentialMasterCell.init();
     //   double u0 = potentialMasterCell.computeAll(true);
      //  System.out.println(kcals.fromSim(u0) + " " + u0);
/*

        if (configFileName != null) {
            ConfigurationFile config = new ConfigurationFile(configFileName);
            config.initializeCoordinates(box);
            BoxImposePbc.imposePBC(box);
        }
        else {
            ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
            configuration.initializeCoordinates(box);
            potentialMasterCell.init();
            double u0 = potentialMasterCell.computeAll(true);
            int numMol = numMoleulesMOP + numMoleulesGas;
            double x = 1;
            ArrayList<Double> energyList = new ArrayList<>();
            int k = 0;
            while (u0 > 1e5 * numMol) {
                // System.out.println(x +" before");
                x *= 0.99;
                // System.out.println( x +" =x");
                for (j = 0; j < pairAtomSize; j++) {
                    p2LJ1[j].setSigmaNew(x * sigmaIJ1[j]);
                    ((P2SoftSphericalSumTruncatedForceShifted) p2lj1[j]).setTruncationRadius(rc);
                }
                u0 = potentialMasterCell.computeAll(false);
                energyList.add(u0);
                k++;
                if (energyList.size() >= 2) {
                    System.out.println(energyList);
                    double diff = energyList.get(k - 1) - energyList.get(k - 2);
                    double percDiff = (diff / (energyList.get(k - 1))) * 100;
                    System.out.println(diff + " " + percDiff);
                    if (percDiff < 1) {
                        break;
                    }
                }
            }

            integratorMC.reset();
            while (u0 > 1e5 * numMol) {
                while (u0 > 1e5 * numMol) {

                    integratorMC.doStep();
                    u0 = integratorMC.getPotentialEnergy();
                    energyList.add(u0);
                    k++;
                    jumperClass(u0, k);
                    // System.out.println("Inside Loop Two - 1");
                }
                while (x < 1 && u0 <= 1e5 * numMol) {
                    //  System.out.println("Inside Loop Two - 2");
                    x /= 0.99;
                    if (x > 1) x = 1;
                    for (j = 0; j < pairAtomSize; j++) {
                        // System.out.println("Inside Loop Two - 3");
                        p2LJ1[j].setSigmaNew(x * sigmaIJ1[j]);
                        ((P2SoftSphericalSumTruncatedForceShifted) p2lj1[j]).setTruncationRadius(rc);
                    }
                    u0 = potentialMasterCell.computeAll(false);
                    //System.out.println(u0 +" inside Array @");
                }
                integratorMC.reset();
            }
        }

*/
    }


    public static void main(String[] args) {
        final String APP_NAME = "Methane Universal";
        GasMOPSimulationParams params = new GasMOPSimulationParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
       /* else {
            params.density = 0.0000003;
            params.configFilename = null; // "octane";
            params.numSteps = 20000;
            params.graphics = true;
        }*/


      //Unit dUnit = new SimpleUnit(Null.DIMENSION, 100, "Density", "g/cm^3", false);
       Unit[] pref = {Dalton.UNIT, Angstrom.UNIT};
        double[] expo = {1.0, -3.0};
        CompoundUnit dalUnit = new CompoundUnit(pref, expo);
        double temperatureK = params.temperatureK;
     //   double temperature = Kelvin.UNIT.toSim(temperatureK);
       // System.out.println("Tsim "+temperature);
        int numMoleculesGas = params.numMoleculesGas;
        int numMoleculesMOP = params.NumMoleculesMOP;
        double density = dalUnit.toSim(params.density);
        long numSteps = params.numSteps;
        String configFilename = params.configFilename;
        double rc = params.rc;
        double pressureKPa = params.pressureKPa;
        Unit pUnit = new PrefixedUnit(Prefix.KILO, Pascal.UNIT);
        double pressure = pUnit.toSim(pressureKPa);
        double rcLimit = params.rcLimit;
        String speciesMOPName = params.speciesMOPName;
        String speciesGasName = params.speciesGasName;
        boolean isDynamic = params.isDynamic;
        boolean dographics = params.graphics;
        boolean isResidual = params.isResidual;
        Vector centreMOP = params.centreMOP;
        Vector boxSize = params.boxSize;

       // System.out.println(numSteps+" steps");
      //  System.out.println("rc: "+rc);
      //  System.out.println("pressure "+ pressureKPa);
       // System.out.println("initial density "+ density);
       // System.out.println("initial density (g/cm^3) "+ density);
        List<Double> rcValues = new ArrayList<>();
        //rcValues.add(rcLimit);
        for(double i=params.tempStart; i<params.tempLimit ; i+=params.tempDiff){
            rcValues.add(i);
        }
        System.out.println(rcValues);
        for (int i=0; i<rcValues.size(); i++){
            double temperature = Kelvin.UNIT.toSim(rcValues.get(i));
            final GasMOPSimulation sim = new GasMOPSimulation(Space3D.getInstance(),speciesMOPName, centreMOP, speciesGasName, density, numMoleculesGas, numMoleculesMOP, temperature, configFilename, rcLimit, pressure, isDynamic, boxSize );
            MeterPotentialEnergyFromIntegrator meterU = new MeterPotentialEnergyFromIntegrator(sim.integratorMC);
            sim.integratorMC.getPotentialCompute().init();
            sim.integratorMC.reset();
            // System.out.println("u0/N "+(meterU.getDataAsScalar()/numMoleculesGas));
            //     Unit kjmol = new UnitRatio(new PrefixedUnit(Prefix.KILO,Joule.UNIT), Mole.UNIT);
            //     System.out.println("u0/N  "+ kjmol.fromSim(meterU.getDataAsScalar() / numMoleculesGas) + " kJ/mol");
            //System.exit(1);

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
            //    System.out.println("Reached before dataforked");
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

            long samples = numSteps / (numMoleculesGas* 8L);
            long bs = samples / 10;
            if (bs == 0) bs = 1;

            MeterWidomInsertion meterWidom = new MeterWidomInsertion(sim.box, sim.random, sim.integratorMC.getPotentialCompute(), sim.integratorMC.getTemperature());
            // meterWidom.setPressure(Bar.UNIT.toSim(15));
            meterWidom.setSpecies(sim.speciesGas);
            meterWidom.setResidual(isResidual);
            meterWidom.setNInsert(100);
            //meterWidom.setPressure(Bar.UNIT.toSim(500));
            meterWidom.setTemperature(sim.integratorMC.getTemperature());
        if(dographics){
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Widom");
            DiameterHashByType dhbt = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
            dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 1);
            //dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 1);
            /*dhbt.setDiameter(sim.speciesMOP.getAtomType(1), 1);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(2), 1);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(3), 1);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 1);*/
            //dhbt.setDiameter(sim.speciesMOP.getAtomType(5), 1);
            //dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 1);
            //dhbt.setDiameter(sim.speciesGrapheneTwo.getAtomType(0), 1);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(0), Color.lightGray);
            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(1), Color.red);
            //   ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(0), Color.cyan);

            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getTypeByName("Cu"), Color.white);
            // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getTypeByName("H"), Color.RED);
            // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getTypeByName("Ar"), Color.cyan);

            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getTypeByName("C_3p"), Color.cyan);
            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getTypeByName("H_p"), Color.cyan);
            // simGraphic.makeAndDisplayFrame("GCMC");
            // ActivityIntegrate ai2 = new ActivityIntegrate(sim.integratorMC);
            // sim.getController().addActivity(ai2, Long.MAX_VALUE, 1.0);
            simGraphic.makeAndDisplayFrame("Widom");
            ActivityIntegrate ai2 = new ActivityIntegrate(sim.integratorMC);
            sim.getController().addActivity(ai2, Long.MAX_VALUE, 1.0);
            return;
        }
            meterWidom.setPositionSource(new RandomPositionSourceRectangular(Space3D.getInstance(), sim.random) {
                public Vector randomPosition() {
                    Vector v;
                    do {
                        v = super.randomPosition();
                    }
                    while (v.getX(0) < 0);
                    return v;
                }
            });
            //   System.out.println("Reached after for loop");
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps/10));
            //sim.integratorMC.getMoveManager().addMCMove(sim.mcMoveMolecule);

            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps/5));
            DataFork muDataFork = new DataFork();
            DataPumpListener muDataPump = new DataPumpListener(meterWidom, muDataFork);
            AccumulatorAverageCollapsing muDataAccu = new AccumulatorAverageCollapsing();
            muDataFork.addDataSink(muDataAccu);
            sim.integratorMC.getEventManager().addListener(muDataPump);
            muDataAccu.setPushInterval(100);
        /*AccumulatorHistory muHistoryA = new AccumulatorHistory(new HistoryCollapsingAverage());
        muDataAccu.addDataSink(muHistoryA);*/
            DataProcessor uProcessor = new DataProcessorFunction(new IFunction() {
                public double f(double x) {
                    if (x == 0) return Double.POSITIVE_INFINITY;
                    return -Math.log(x) * sim.integratorMC.getTemperature();
                }
            });
            muDataFork.addDataSink(uProcessor);

            //  System.out.println("Reached after interval");

/*
       AccumulatorAverageFixed accU = new AccumulatorAverageFixed((numSteps/10)/1000);
        DataPumpListener pumpU = new DataPumpListener(meterU, accU, 500);
        sim.integratorMC.getEventManager().addListener(pumpU);

        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(bs);
        forkP.addDataSink(accP);
        AccumulatorAverageFixed accZ = new AccumulatorAverageFixed(bs);
        dpZ.addDataSink(accZ);
        AccumulatorAverageFixed accZm1oR = new AccumulatorAverageFixed(bs);
        dpZm1oR.addDataSink(accZm1oR);
        DataPumpListener pumpP = new DataPumpListener(meterP, forkP, 8*numMoleculesGas);
        sim.integratorMC.getEventManager().addListener(pumpP);
*/
            long t1 = System.nanoTime();
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps));
            long t2 = System.nanoTime();
            // System.out.println("Done with time");
            IData dataWidom = muDataAccu.getData();
            double scalar = meterWidom.getDataAsScalar();
            double avg = dataWidom.getValue(muDataAccu.AVERAGE.index);
            double errRho = muDataAccu.getData(muDataAccu.ERROR).getValue(0);
            double corRho = muDataAccu.getData(muDataAccu.BLOCK_CORRELATION).getValue(0);
            double chemicalPot = -Constants.BOLTZMANN_K*temperatureK*Math.log(avg);
            double chemPot10 = -Constants.BOLTZMANN_K*temperatureK*Math.log10(avg);
            System.out.println("conf Names " +params.speciesMOPName + " " + params.speciesGasName);
            System.out.println("temp "+ rcValues.get(i)+ " boxsize "+ boxSize  );
            System.out.println(avg + " widom " + chemicalPot + " "+ chemPot10 + " " + errRho +  " " +corRho);
            double[] arrBox = boxSize.toArray();
            System.out.println("B2 " + -0.5*arrBox[0]*arrBox[1]*arrBox[2]*(1-avg) + " " + -0.5*arrBox[0]*arrBox[1]*arrBox[2]*(avg-1));

       /* IData dataU = muDataAccu.getData();
        double avgU = dataU.getValue(AccumulatorAverage.AVERAGE.index) / numMoleculesGas;
        double errU = dataU.getValue(AccumulatorAverage.ERROR.index) / numMoleculesGas;
        double corU = dataU.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("U: " +" "+avgU+ "   err: "+" "+errU+"   cor: " +" "+corU);

        UnitRatio jouleMole = new UnitRatio(Joule.UNIT, Mole.UNIT);
        double valjouleMole = jouleMole.fromSim(dataU.getValue(AccumulatorAverage.AVERAGE.index/numMoleculesGas));
        System.out.println(valjouleMole + " J/mol");*/
     /*   IData dataP = accP.getData();
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
        System.out.println("(Z-1)/rho: "+" "+avgZ_+"   err: "+" "+errZ_+"   cor: "+" "+corZ_);*/

            System.out.println("time: "+" "+(t2-t1)/1e9);
        }
    }

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



    public static class GasMOPSimulationParams extends ParameterBase {
        public String speciesMOPName = "F://Avagadro//molecule//co2";
        public String speciesGasName = "F://Avagadro//molecule//Ar" ;
        public double temperatureK = 330;
        public int numMoleculesGas = 1;
        public int NumMoleculesMOP = 1;
        //public int pressure = 10;
        public double density = 1e-9;
        public boolean graphics = true;
        public long numSteps = 1000000;
        public Vector boxSize = new Vector3D(14,14,14);
        public String configFilename = null;
        public double rc = 10;
       // public double rcLimit = boxSize.toArray()[0]/2 + 1 ;
        public double rcLimit =10 ;
        public double stepsize = 3;
        public Vector centreMOP = new Vector3D(0.0,0.0, 0.0);

        public double tempStart = 250;
        public double tempDiff = 20;
        public double tempLimit = 375;
        public double pressureKPa = 1402;
        public boolean isDynamic = true;
        public boolean isResidual = true;

    }
}
