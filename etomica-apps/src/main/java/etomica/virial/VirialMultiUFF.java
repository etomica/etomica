package etomica.virial;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.potential.UFF.*;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.units.dimensions.*;
import etomica.units.dimensions.Dimension;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.collections.IntArrayList;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.*;
import etomica.virial.simulations.SimulationVirialOverlap2;
import etomica.virial.wheatley.ClusterWheatleySoftMix;


import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;

public class VirialMultiUFF {
    public static ISpecies species1, species2;
    public static Box box;
    public static int atom1, atom2, atom3, i=0;
    public static String atomName1, atomName2, atomName3;


    public static void main(String[] args) {
        int i;
        VirialMixedSCParam params = new VirialMixedSCParam();
        final int[] nTypes = params.nTypes;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double refFreq = params.refFrac;
        final int nPoints = params.nTypes[0] + params.nTypes[1];
        double sigmaHSRef = params.sigmaHSRef;
        String confName1 = params.confName1;
        String confName2 = params.confName2;
        Space space = Space3D.getInstance();
        boolean isCommandline = args.length > 0;

        PDBReader pdbReader1 = new PDBReader();
        species1 = PDBReader.getSpecies(confName1);
        System.out.println(species1);
        printCOM(species1);
        //List<AtomType> atomTypes1 = species1.getUniqueAtomTypes();
        List<List<AtomType>> pairsAtoms1 = getSpeciesPairs(species1);


        PDBReaderReplica pdbReader2 = new PDBReaderReplica();
        int pairAtomSize = pairsAtoms1.size();

        species2 = PDBReaderReplica.getSpecies(confName2);
        printCOM(species2);
        List<List<AtomType>> pairsAtoms2 = getSpeciesPairs(species2);
        int pairAtomSize2 = pairsAtoms2.size();
        SpeciesManager sm1 = new SpeciesManager.Builder().addSpecies(species1).addSpecies(species2).build();
        UniversalSimulation.makeAtomPotentials(sm1);
        System.out.println(nPoints+" at "+temperature+"K");
        temperature = Kelvin.UNIT.toSim(temperature);

        ArrayList<ArrayList<Integer>> connectedAtoms1 = PDBReader.getConnectivityWithoutRunning();
        System.out.println(connectedAtoms1);
        ArrayList<ArrayList<Integer>> connectivityModified1 = PDBReader.getConnectivityModifiedWithoutRunning();
        System.out.println(connectivityModified1);
        Map<Integer,String> atomMap1 = PDBReader.getAtomMapWithoutRunning();
        System.out.println(atomMap1);
        HashMap<Integer, String> atomMapModified1 = PDBReader.getAtomMapModifiedWithoutRunning();
        System.out.println(atomMapModified1);
        ArrayList<Integer> bondList1 = PDBReader.getBondList(connectedAtoms1, atomMap1);
        System.out.println(bondList1);
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        Map<String, double[]> atomicPotMap1 = PDBReader.atomicPotMap();
        System.out.println(atomicPotMap1);
        ArrayList<Integer> bondsNum1 = PDBReader.getBonds();

        //Map<Integer, String> atomIdentifierMapModified1 = PDBReader.atomIdentifierMapModified(connectivityModified1, atomMapModified1);
        Map<Integer, String> atomIdentifierMapModified1 = PDBReader.getModifiedAtomIdentifierMap();
        List<int[]>dupletsSorted1= PDBReader.getDupletesSorted();
        List<int[]>tripletsSorted1= PDBReader.getAnglesSorted();
        List<int[]>quadrupletsSorted1= PDBReader.getTorsionSorted();

        Map<String[],List<int[]>> bondTypesMap1= PDBReader.idenBondTypes(dupletsSorted1, atomIdentifierMapModified1);
        Map<String[],List<int[]>> angleTypesMap1= PDBReader.idenAngleTypes(tripletsSorted1, atomIdentifierMapModified1);
        Map<String[],List<int[]>> torsionTypesMap1= PDBReader.idenTorsionTypes(quadrupletsSorted1, atomIdentifierMapModified1);
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
            for (int j = 0; j < innerList.size(); j++) {
                intArrayList.add(innerList.get(j));
                System.out.println(intArrayList);
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
        //System.out.println(modifiedOutput1);
        //System.out.println(Arrays.toString(dupletsIntArrayList1));

        box = new Box(space);
        // FIll it up with Meyer Function and Target Diagram
        PotentialMoleculePair pTargetA = new PotentialMoleculePair(space, sm1);
        PotentialMoleculePair pTargetB = new PotentialMoleculePair(space, sm1);
        PotentialMoleculePair pTargetAB = new PotentialMoleculePair(space, sm1);
        MayerGeneral fTargetA = new MayerGeneral(pTargetA);
        MayerGeneral fTargetB = new MayerGeneral(pTargetB);
        MayerGeneral fTargetAB = new MayerGeneral(pTargetAB);
        MayerFunction[][] allF = new MayerFunction[][]{{fTargetA,fTargetAB},{fTargetAB,fTargetB}};
        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }
            @Override
            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef * sigmaHSRef ? 1 : 0;
            }
        };
        boolean alkaneFlex = nPoints > 2;
        boolean[] boolOne = {false, false};
        boolean[] boolTwo = {true, true};
        VirialDiagramsMix alkaneDiagrams = new VirialDiagramsMix(nPoints, boolOne, boolTwo);
        alkaneDiagrams.setDoReeHoover(false);
        System.out.println(allF.length);
        ClusterAbstract targetCluster = new ClusterWheatleySoftMix(nPoints, nTypes, allF, 1e-12);
        ClusterChainHS refCluster = new ClusterChainHS(nPoints, fRefPos);

        double vhs = (4.0 / 3.0) * Math.PI * params.sigmaHSRef * params.sigmaHSRef *params.sigmaHSRef;
        final double refIntegral = SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1);


        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
        int[] targetDiagramNumbers = new int[0];
        boolean[] diagramFlexCorrection = null;

        targetDiagramNumbers = new int[targetDiagrams.length];
        System.out.println("individual clusters:");
        diagramFlexCorrection = new boolean[targetDiagrams.length];

        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        for (i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }

        System.out.println("sigmaHSRef: "+params.sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+refIntegral);

        PotentialMasterBonding.FullBondingInfo bondingInfo1 = new PotentialMasterBonding.FullBondingInfo(sm1);
        doBondStrech(species1, bondTypesMap1, angleTypesMap1, torsionTypesMap1,bondsNum1,bondList1, quadrupletsSorted1, atomIdentifierMapModified1,atomicPotMap1, bondingInfo1);
        LJUFF[] p2LJ = new LJUFF[pairAtomSize];
        doLJ(pairsAtoms1,pTargetA, p2LJ);

        temperature = Kelvin.UNIT.toSim(temperature);
        ArrayList<ArrayList<Integer>> connectedAtoms2 =PDBReaderReplica.getConnectivityWithoutRunning();
        ArrayList<ArrayList<Integer>> connectivityModified2 = PDBReaderReplica.getConnectivityModifiedWithoutRunning();
        Map<Integer,String> atomMap2 = PDBReaderReplica.getAtomMapWithoutRunning();
        HashMap<Integer, String> atomMapModified2 = PDBReaderReplica.getAtomMapModifiedWithoutRunning();
        ArrayList<Integer> bondList2 = PDBReaderReplica.getBondList(connectedAtoms2, atomMap2);
       // System.out.println(bondList2 +" bondList1");
        Map<String, double[]> atomicPotMap2 = PDBReaderReplica.atomicPotMap();
        Map<Integer, String> atomIdentifierMapModified2 = PDBReaderReplica.getatomIdentifierMapModified();
       // System.out.println(atomIdentifierMapModified2 + "atomIdentifierMapModified1");
        List<int[]>dupletsSorted2= PDBReaderReplica.getDupletesSorted();
        List<int[]>tripletsSorted2=PDBReaderReplica.getAnglesSorted();
        List<int[]>quadrupletsSorted2=PDBReaderReplica.getTorsionSorted();
        ArrayList<Integer> bondsNum2 = PDBReaderReplica.getBonds();
       // System.out.println(Arrays.deepToString(dupletsSorted2.toArray()) + " dupletssorted");
       // System.out.println(atomIdentifierMapModified2);

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
            for (int j = 0; j < innerList.size(); j++) {
                intArrayList.add(innerList.get(j));
                System.out.println(intArrayList);
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
        //System.out.println(modifiedOutput2);
        //System.out.println(Arrays.toString(dupletsIntArrayList2));
        //targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        for (i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }
        System.out.println("sigmaHSRef: "+params.sigmaHSRef);
        // eovererr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);
        doBondStrech(species2,bondTypesMap2, angleTypesMap2, torsionTypesMap2, bondsNum2, bondList2,quadrupletsSorted2, atomIdentifierMapModified2, atomicPotMap2, bondingInfo1);
        LJUFF[] p2LJ2 = new LJUFF[pairAtomSize2];
        doLJ(pairsAtoms2, pTargetB,  p2LJ2);

        System.out.println(steps+" steps (100 blocks of "+steps/1000+")");
        System.out.println(species1);
        System.out.println(species2);
        List<AtomType> list1 = species1.getUniqueAtomTypes();
        List<AtomType> list2 = species2.getUniqueAtomTypes();
        List<AtomType> list3 = new ArrayList<>(list1);
        int list2Size = list2.size();
        boolean isEqual =false;
        for(i=0; i<list2Size; i++) {
            String name = list2.get(i).getName();
            for(int j =0; j<list3.size(); j++){
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
            for (int j = 0; j < list3.size(); j++) {
                if(i<=j){
                    List<AtomType> subPair = new ArrayList<>();
                    subPair.add(list3.get(i));
                    subPair.add(list3.get(j));
                    pairsAtomsTotal.add(subPair);
                }
            }
        }

        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species1,species2}, nTypes, temperature, refCluster, targetCluster);
        sim.setExtraTargetClusters(targetDiagrams);
        sim.setBondingInfo(bondingInfo1);
        sim.setIntraPairPotentials(pTargetA.getAtomPotentials());
        sim.setIntraPairPotentials(pTargetB.getAtomPotentials());
        sim.setRandom(new RandomMersenneTwister(2));
        sim.init();
        int[] seeds = sim.getRandomSeeds();
        System.out.println("Random seeds: "+ Arrays.toString(seeds));
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        PotentialCompute pc0 = sim.integrators[0].getPotentialCompute();
        MCMoveClusterStretch mcMoveStretchRef1 = new MCMoveClusterStretch(pc0, space, dupletsIntArrayList1, sim.getRandom(), 0.05);
        mcMoveStretchRef1.setSpecies(species1);
        MCMoveClusterAngle mcMoveBendRef1 = new MCMoveClusterAngle(pc0, space, dupletsIntArrayList1, sim.getRandom(), 0.05);
        mcMoveBendRef1.setSpecies(species1);
        sim.integrators[0].getMoveManager().addMCMove(mcMoveStretchRef1);
        sim.integrators[0].getMoveManager().addMCMove(mcMoveBendRef1);
        MCMoveClusterStretch mcMoveStretchRef2 = new MCMoveClusterStretch(pc0, space, dupletsIntArrayList2, sim.getRandom(), 0.05);
        mcMoveStretchRef2.setSpecies(species2);
        MCMoveClusterAngleBend mcMoveBendRef2 = new MCMoveClusterAngleBend(pc0, sim.getRandom(), 0.05, space);
        mcMoveBendRef2.setSpecies(species2);
        sim.integrators[0].getMoveManager().addMCMove(mcMoveStretchRef2);
        sim.integrators[0].getMoveManager().addMCMove(mcMoveBendRef2);
        PotentialCompute pc1 = sim.integrators[1].getPotentialCompute();
        MCMoveClusterStretch mcMoveStretchTarget1 = new MCMoveClusterStretch(pc1, space, dupletsIntArrayList1, sim.getRandom(), 0.05);
        mcMoveStretchTarget1.setSpecies(species1);
        MCMoveClusterAngleBend mcMoveBendTarget1 = new MCMoveClusterAngleBend(pc1, sim.getRandom(), 0.05, space);
        mcMoveBendTarget1.setSpecies(species1);
        sim.integrators[1].getMoveManager().addMCMove(mcMoveStretchTarget1);
        sim.integrators[1].getMoveManager().addMCMove(mcMoveBendTarget1);
        MCMoveClusterStretch mcMoveStretchTarget2 = new MCMoveClusterStretch(pc1, space, dupletsIntArrayList2, sim.getRandom(), 0.05);
        mcMoveStretchTarget2.setSpecies(species2);
        MCMoveClusterAngleBend mcMoveBendTarget2 = new MCMoveClusterAngleBend(pc1, sim.getRandom(), 0.05, space);
        mcMoveBendTarget2.setSpecies(species2);
        sim.integrators[1].getMoveManager().addMCMove(mcMoveStretchTarget2);
        sim.integrators[1].getMoveManager().addMCMove(mcMoveBendTarget2);
        if (alkaneFlex) {
            int[] constraintMap = new int[nPoints+1];
            for (i=0; i<nPoints; i++) {
                constraintMap[i] = i;
            }
            constraintMap[nPoints] = 0;
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[0]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).setConstraintMap(constraintMap);
        }

        if (refFreq >= 0) {
            sim.integratorOS.setAdjustStepFraction(false);
            sim.integratorOS.setRefStepFraction(refFreq);
        }

        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;

        int pairAtom1Size = pairsAtoms1.size();
        int paitAtom2Size = pairsAtoms2.size();
        int pairAtomsTotalSize = pairsAtomsTotal.size();
        LJUFF[] p2LJAB = new LJUFF[pairAtomsTotalSize];
        doLJ(pairsAtomsTotal, pTargetAB, p2LJAB);


        if(false) {
            double size =20;
            sim.box[0].getBoundary().setBoxSize(etomica.space.Vector.of(new double[]{size, size, size}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox1.setPixelUnit(new Pixel(300.0 / size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys) displayBox1.canvas).setBackgroundColor(Color.WHITE);


            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box[0], sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box[1], sim.getRandom());
            displayBox1.setColorScheme(colorScheme);


            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 10, false);
            sim.equilibrate(null, 20, false);
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }

            final DisplayTextBox averageBox = new DisplayTextBox();
            averageBox.setLabel("Average");
            final DisplayTextBox errorBox = new DisplayTextBox();
            errorBox.setLabel("Error");
            JLabel jLabelPanelParentGroup = new JLabel("B" + nPoints + " (L/mol)^" + (nPoints - 1));
            final JPanel panelParentGroup = new JPanel(new BorderLayout());
            panelParentGroup.add(jLabelPanelParentGroup, Constants.CompassDirection.NORTH.toString());
            panelParentGroup.add(averageBox.graphic(), BorderLayout.WEST);
            panelParentGroup.add(errorBox.graphic(), BorderLayout.EAST);
            simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());

            IAction pushAnswer = new IAction() {
                DataDouble data = new DataDouble();

                public void actionPerformed() {
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    data.x = ratio;
                    averageBox.putData(data);
                    data.x = error;
                    errorBox.putData(data);
                }
            };
            IDataInfo dataInfo = new DataDouble.DataInfoDouble("B" + nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints - 1}));
            Unit unit = new CompoundUnit(new Unit[]{new UnitRatio(Liter.UNIT, Mole.UNIT)}, new double[]{nPoints - 1});
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            averageBox.setUnit(unit);
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            errorBox.setUnit(unit);
            sim.integratorOS.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));
            return;
        }

        long t1 = System.nanoTime();
        // if running interactively, don't use the file
        String refFileName = isCommandline ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        System.out.println("Hello");
        sim.initRefPref(refFileName, steps/40);


        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, 1000);
        sim.setAccumulatorBlockSize(steps);

        System.out.println("equilibration finished");
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int)steps);
        System.out.println("MC Move step sizes (ref)    "
                +(sim.mcMoveRotate==null ? "" : (""+sim.mcMoveRotate[0].getStepSize()))+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +(sim.mcMoveRotate==null ? "" : (""+sim.mcMoveRotate[1].getStepSize()))+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[1].getStepSize())));

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(ai);
        long t2 = System.nanoTime();



        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        String[] extraNames = new String[targetDiagrams.length];
        for (i=0; i<targetDiagrams.length; i++) {
            String n = "";
            if (targetDiagramNumbers[i] < 0) {
                n = "diagram " + (-targetDiagramNumbers[i]) + "bc";
            } else {
                n = "diagram " + targetDiagramNumbers[i];

                if (diagramFlexCorrection[i]) {
                    n += "c";
                }
            }
            extraNames[i] = n;
        }
        sim.printResults(refIntegral, extraNames);
        System.out.println("time: "+(t2-t1)/1e9);

    }
    public static void doLJ(List<List<AtomType>> pairsAtoms, PotentialMoleculePair pTarget, LJUFF[] p2LJ){
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
            //sigmaIJ[i] = (sigmaIKey + sigmaJKey) / 2;
            p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);
            System.out.println(atomTypeOne + " " + atomTypeTwo +" "+ Arrays.toString(iKey) +" " +Arrays.toString(jKey));
            pTarget.setAtomPotentials(atomNameOne, atomNameTwo, p2LJ[i], new double[]{1, 0, 0, 1});
            i++;
        }
    }
    public static void doNoble(ISpecies species1, Map<Integer, String> atomIdentifierMapModified1,Map<String, double[]> atomicPotMap1, PotentialMasterBonding.FullBondingInfo bondingInfo1 ){

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
                bondingInfo1.setBondingPotentialPair(species1, p2Bond, bonds );
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
    public static List<List<AtomType>> getSpeciesPairs (ISpecies species){
        List<List<AtomType>> pairsAtoms1 = new ArrayList<>();
        List<AtomType> atomTypes1 = species.getUniqueAtomTypes();
        int i, j;
        for(i=0; i<atomTypes1.size(); i++) {
            for (j = 0; j < atomTypes1.size(); j++) {
                if(i<=j){
                    List<AtomType> subPair = new ArrayList<>();
                    subPair.add(species.getAtomType(i));
                    subPair.add(species.getAtomType(j));
                    pairsAtoms1.add(subPair);
                }
            }
        }
        return pairsAtoms1;
    }
    public static void printCOM(ISpecies species){
        Space space = Space3D.getInstance();
        Vector center = space.makeVector();
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
       // System.out.println(center + " 2 out");
    }

    public static class VirialMixedSCParam extends ParameterBase {
        String confName1 = "F://Avagadro//molecule//h2";
        String confName2 = "F://Avagadro//molecule//o2";
        // don't change these
        public int[] nTypes = new int[]{1, 1};
        public double temperature = 300;
        public long numSteps =100000;
        public double refFrac = -1;
        public double sigmaHSRef = 5.23;
        //public Level level = Level.CLASSICAL;
        public boolean doHist = false;
      //  public Nonadditive nonAdditive = Nonadditive.NONE;
        public boolean useSZ = false;
        public int rc = 10;
    }
}
