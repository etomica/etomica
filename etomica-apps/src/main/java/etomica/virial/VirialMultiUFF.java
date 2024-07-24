package etomica.virial;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.potential.UFF.*;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.util.ParameterBase;
import etomica.util.collections.IntArrayList;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.*;
import etomica.virial.simulations.SimulationVirialOverlap2;
import etomica.virial.wheatley.ClusterWheatleySoftMix;
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
       // double temperature = params.temperature;
        //long steps = params.numSteps;
        double refFreq = params.refFrac;
        final int nPoints = params.nTypes[0] + params.nTypes[1];
        double sigmaHSRef = params.sigmaHSRef;
        //String confName1 = params.confName1;
       // String confName2 = params.confName2;
        Vector centreMol = params.centreMoleculeOne;
        Space space = Space3D.getInstance();
        boolean isCommandline = args.length > 0;
        boolean doTrunc = params.doTrunc;
        List<Double> truncRad = new ArrayList<>();
        List<String> moleName = new ArrayList<>();
        moleName.add(params.confListName1);
        moleName.add(params.confListName2);
        moleName.add(params.confListName3);
        moleName.add(params.confListName5);
        moleName.add(params.confListName4);
        List<String> moleName2 = new ArrayList<>();
        moleName2.add(params.confName1);
     //   moleName2.add(params.confName2);
     //   moleName2.add(params.confName3);
      //  moleName2.add(params.confName4);

      /*  for ( double j= params.start; j<params.truncLimit; j+=params.truncDiff ){
            truncRad.add(j);
        }*/
      /*  for ( double j= params.tempStart; j<params.tempFinal; j+=params.tempDiff ){
            truncRad.add(j);
        }*/
       // System.out.println(truncRad);
        double rc = params.start;
        System.out.println(rc);
      /*  for (int j=0; j< moleName.size(); j++){
            System.out.println(moleName.get(j));
        }*/
    //    System.exit(1);
        for(int p=0; p< moleName.size(); p++){
        for ( int m=0; m< moleName2.size(); m++){
        //    double rc = truncRad.get(m);
            double temperatureSim = params.temperature;
            PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
            String confName2 = moleName.get(p);
            String confName1 = moleName2.get(m);
            species1 = pdbReaderMOP.getSpecies(confName1, false, centreMol, false);
            //System.out.println(species1);
            printCOM(species1);
            //List<AtomType> atomTypes1 = species1.getUniqueAtomTypes();
            List<List<AtomType>> pairsAtoms1 = getSpeciesPairs(species1);
            int pairAtomSize = pairsAtoms1.size();

            PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
            species2 = pdbReaderReplica.getSpecies(confName2, true, new Vector3D(0,0,0), false);
           // printCOM(species2);
            List<List<AtomType>> pairsAtoms2 = getSpeciesPairs(species2);
            int pairAtomSize2 = pairsAtoms2.size();
            SpeciesManager sm1 = new SpeciesManager.Builder().addSpecies(species1).addSpecies(species2).build();
            UniversalSimulation.makeAtomPotentials(sm1);
            //System.out.println(nPoints+" at "+temperature+"K");
            //temperature = Kelvin.UNIT.toSim(temperature);

            ArrayList<ArrayList<Integer>> connectedAtoms1 = pdbReaderMOP.getConnectivity();
            // System.out.println(connectedAtoms1);
            ArrayList<ArrayList<Integer>> connectivityModified1 = pdbReaderMOP.getConnectivityModified();
            //System.out.println(connectivityModified1);
            Map<Integer,String> atomMap1 = pdbReaderMOP.getAtomMap();
            // System.out.println(atomMap1);
            Map<Integer, String> atomMapModified1 = pdbReaderMOP.getAtomMapModifiedWithoutRunning();
            //System.out.println(atomMapModified1);
            ArrayList<Integer> bondList1 = pdbReaderMOP.getBondList(connectedAtoms1, atomMap1);
            // System.out.println(bondList1);
            Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
            Map<String, double[]> atomicPotMap1 = pdbReaderMOP.atomicPotMap();
            //System.out.println(atomicPotMap1);
            ArrayList<Integer> bondsNum1 = pdbReaderMOP.getBonds();

            //Map<Integer, String> atomIdentifierMapModified1 = PDBReader.atomIdentifierMapModified(connectivityModified1, atomMapModified1);
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
            // System.out.println(modifiedOutput1 +" modified");
            IntArrayList[] dupletsIntArrayList1 = new IntArrayList[modifiedOutput1.size()];

            for (i = 0; i < modifiedOutput1.size(); i++) {
                ArrayList<Integer> innerList = modifiedOutput1.get(i);
                IntArrayList intArrayList = new IntArrayList(innerList.size());
                for (int j = 0; j < innerList.size(); j++) {
                    intArrayList.add(innerList.get(j));
                    // System.out.println(intArrayList);
                }
                dupletsIntArrayList1[i] = intArrayList;
                //  System.out.println(dupletsIntArrayList[i]);
            }
            for (IntArrayList list : dupletsIntArrayList1) {
                for (i = 0; i < list.size(); i++) {
                    int value = list.getInt(i);
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
            // System.out.println(allF.length);
            ClusterAbstract targetCluster = new ClusterWheatleySoftMix(nPoints, nTypes, allF, 1e-12);
            ClusterChainHS refCluster = new ClusterChainHS(nPoints, fRefPos);

            double vhs = (4.0 / 3.0) * Math.PI * params.sigmaHSRef * params.sigmaHSRef *params.sigmaHSRef;
            final double refIntegral = SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1);


            ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
            int[] targetDiagramNumbers = new int[0];
            boolean[] diagramFlexCorrection = null;

            targetDiagramNumbers = new int[targetDiagrams.length];
            //System.out.println("individual clusters:");
            diagramFlexCorrection = new boolean[targetDiagrams.length];



         //   System.out.println("sigmaHSRef: "+params.sigmaHSRef);
          //  System.out.println("B"+nPoints+"HS: "+refIntegral);

            PotentialMasterBonding.FullBondingInfo bondingInfo1 = new PotentialMasterBonding.FullBondingInfo(sm1);
            doBondStrech(species1, bondTypesMap1, angleTypesMap1, torsionTypesMap1,bondsNum1,bondList1, quadrupletsSorted1, atomIdentifierMapModified1,atomicPotMap1, bondingInfo1);
            LJUFF[] p2LJ = new LJUFF[pairAtomSize];
            IPotential2[] p2lj = new IPotential2[pairAtomSize];
            doLJ(pairsAtoms1,pTargetA, p2LJ, p2lj,rc, doTrunc);

            // temperature = Kelvin.UNIT.toSim(temperature);
            ArrayList<ArrayList<Integer>> connectedAtoms2 =pdbReaderReplica.getConnectivityWithoutRunning();
            ArrayList<ArrayList<Integer>> connectivityModified2 = pdbReaderReplica.getConnectivityModifiedWithoutRunning();
            Map<Integer,String> atomMap2 = pdbReaderReplica.getAtomMapWithoutRunning();
            HashMap<Integer, String> atomMapModified2 = pdbReaderReplica.getAtomMapModifiedWithoutRunning();
            ArrayList<Integer> bondList2 = pdbReaderReplica.getBondList(connectedAtoms2, atomMap2);
            // System.out.println(bondList2 +" bondList1");
            Map<String, double[]> atomicPotMap2 = pdbReaderReplica.atomicPotMap();
            Map<Integer, String> atomIdentifierMapModified2 = pdbReaderReplica.getatomIdentifierMapModified();
            // System.out.println(atomIdentifierMapModified2 + "atomIdentifierMapModified1");
            List<int[]>dupletsSorted2= pdbReaderReplica.getDupletesSorted();
            List<int[]>tripletsSorted2=pdbReaderReplica.getAnglesSorted();
            List<int[]>quadrupletsSorted2=pdbReaderReplica.getTorsionSorted();
            ArrayList<Integer> bondsNum2 = pdbReaderReplica.getBonds();
            // System.out.println(Arrays.deepToString(dupletsSorted2.toArray()) + " dupletssorted");
            // System.out.println(atomIdentifierMapModified2);

            Map<String[],List<int[]>> bondTypesMap2= pdbReaderReplica.idenBondTypes(dupletsSorted2, atomIdentifierMapModified2);
            Map<String[],List<int[]>> angleTypesMap2= pdbReaderReplica.idenAngleTypes(tripletsSorted2, atomIdentifierMapModified2);
            Map<String[],List<int[]>> torsionTypesMap2= pdbReaderReplica.idenTorsionTypes(quadrupletsSorted2, atomIdentifierMapModified2);
            // System.out.println(connectivityModified2);

            ArrayList<ArrayList<Integer>> modifiedOutput2 = new ArrayList<>();
            for (ArrayList<Integer> innerList : connectivityModified2) {
                ArrayList<Integer> modifiedInnerList = new ArrayList<>(innerList.subList(1, innerList.size()));
                modifiedOutput2.add(modifiedInnerList);
            }
            // System.out.println(modifiedOutput2 +" modified");
            IntArrayList[] dupletsIntArrayList2 = new IntArrayList[modifiedOutput2.size()];


            for (i = 0; i < modifiedOutput2.size(); i++) {
                ArrayList<Integer> innerList = modifiedOutput2.get(i);
                IntArrayList intArrayList = new IntArrayList(innerList.size());
                for (int j = 0; j < innerList.size(); j++) {
                    intArrayList.add(innerList.get(j));
                    //  System.out.println(intArrayList);
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



            doBondStrech(species2,bondTypesMap2, angleTypesMap2, torsionTypesMap2, bondsNum2, bondList2,quadrupletsSorted2, atomIdentifierMapModified2, atomicPotMap2, bondingInfo1);
            LJUFF[] p2LJ2 = new LJUFF[pairAtomSize2];
            IPotential2[] p2lj2 = new IPotential2[pairAtomSize2];
            doLJ(pairsAtoms2, pTargetB,  p2LJ2, p2lj2, rc, doTrunc);

            // System.out.println(steps+" steps (100 blocks of "+steps/1000+")");
            //System.out.println(species1);
            // System.out.println(species2);
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
            //    System.out.println(list3);
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
            int pairAtom1Size = pairsAtoms1.size();
            int paitAtom2Size = pairsAtoms2.size();
            int pairAtomsTotalSize = pairsAtomsTotal.size();
            LJUFF[] p2LJAB = new LJUFF[pairAtomsTotalSize];
            IPotential2[] p2ljAB = new IPotential2[pairAtomsTotalSize];
            doLJ(pairsAtomsTotal, pTargetAB, p2LJAB, p2ljAB, rc, doTrunc);

            // System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");
            long t1Main = System.nanoTime();
            System.out.println();
            //  for( int j=1; j<11;j++ ){
            double temperature = Kelvin.UNIT.toSim(temperatureSim);
            long steps = params.numSteps;
            final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species1,species2}, nTypes, temperature, refCluster, targetCluster);
            targetCluster.setTemperature(temperature);
            refCluster.setTemperature(temperature);
            for (i=0; i<targetDiagrams.length; i++) {
                targetDiagrams[i].setTemperature(temperature);
            }
            refCluster.setTemperature(temperature);
            for (i=0; i<targetDiagrams.length; i++) {
                targetDiagrams[i].setTemperature(temperature);
            }
            sim.setExtraTargetClusters(targetDiagrams);
            sim.setBondingInfo(bondingInfo1);
            sim.setIntraPairPotentials(pTargetA.getAtomPotentials());
            sim.setIntraPairPotentials(pTargetB.getAtomPotentials());
            sim.setRandom(new RandomMersenneTwister(2));
            sim.init();
            int[] seeds = sim.getRandomSeeds();
            //System.out.println("Random seeds: "+ Arrays.toString(seeds));
            sim.integratorOS.setAggressiveAdjustStepFraction(true);
            PotentialCompute pc0 = sim.integrators[0].getPotentialCompute();
            MCMoveClusterStretch mcMoveStretchRef1 = new MCMoveClusterStretch(pc0, space, dupletsIntArrayList1, sim.getRandom(), 0.05);
            mcMoveStretchRef1.setSpecies(species1);
            MCMoveClusterAngleBend mcMoveBendRef1 = new MCMoveClusterAngleBend(pc0, sim.getRandom(), 0.05, space);
            mcMoveBendRef1.setSpecies(species1);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveStretchRef1);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveBendRef1);
            // MCMoveClusterMoleculeMulti mcMoveClusterMultiRef2 = new MCMoveClusterMoleculeMulti(sim.getRandom(), box);
            // sim.integrators[0].getMoveManager().addMCMove(mcMoveClusterMultiRef2);
            //MCMoveClusterStretch mcMoveStretchRef2 = new MCMoveClusterStretch(pc0, space, dupletsIntArrayList2, sim.getRandom(), 0.05);
            //mcMoveStretchRef2.setSpecies(species2);
            //MCMoveClusterAngleBend mcMoveBendRef2 = new MCMoveClusterAngleBend(pc0, sim.getRandom(), 0.05, space);
            //mcMoveBendRef2.setSpecies(species2);
            //sim.integrators[0].getMoveManager().addMCMove(mcMoveStretchRef2);
            //sim.integrators[0].getMoveManager().addMCMove(mcMoveBendRef2);
            PotentialCompute pc1 = sim.integrators[1].getPotentialCompute();
            MCMoveClusterStretch mcMoveStretchTarget1 = new MCMoveClusterStretch(pc1, space, dupletsIntArrayList1, sim.getRandom(), 0.05);
            mcMoveStretchTarget1.setSpecies(species1);
            MCMoveClusterAngleBend mcMoveBendTarget1 = new MCMoveClusterAngleBend(pc1, sim.getRandom(), 0.05, space);
            mcMoveBendTarget1.setSpecies(species1);
            sim.integrators[1].getMoveManager().addMCMove(mcMoveStretchTarget1);
            sim.integrators[1].getMoveManager().addMCMove(mcMoveBendTarget1);
            //MCMoveClusterMoleculeMulti mcMoveClusterMultiTarget2 = new MCMoveClusterMoleculeMulti(sim.getRandom(), box);
            //sim.integrators[1].getMoveManager().addMCMove(mcMoveClusterMultiTarget2);
            //MCMoveClusterStretch mcMoveStretchTarget2 = new MCMoveClusterStretch(pc1, space, dupletsIntArrayList2, sim.getRandom(), 0.05);
            //mcMoveStretchTarget2.setSpecies(species2);
            //MCMoveClusterAngleBend mcMoveBendTarget2 = new MCMoveClusterAngleBend(pc1, sim.getRandom(), 0.05, space);
            //mcMoveBendTarget2.setSpecies(species2);
            //sim.integrators[1].getMoveManager().addMCMove(mcMoveStretchTarget2);
            //sim.integrators[1].getMoveManager().addMCMove(mcMoveBendTarget2);
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
            System.out.println(confName1 +" " +confName2 + " " +rc);
            sim.integratorOS.setNumSubSteps(1000);
            sim.integratorOS.setAggressiveAdjustStepFraction(true);
          //  System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
            //System.out.println(doTrunc);
            steps /= 1000;




       /* if(false) {
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
        }*/

            long t1 = System.nanoTime();
            // if running interactively, don't use the file
            // String refFileName = isCommandline ? "refpref"+nPoints+"_"+temperature : null;
            // this will either read the refpref in from a file or run a short simulation to find it
            //System.out.println("Hello");
            // sim.initRefPref(refFileName, steps/40);


            // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
            // if it does continue looking for a pref, it will write the value to the file
            // sim.equilibrate(refFileName, steps/20);
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
            // System.out.println("time: "+(t2-t1)/1e9);
            System.out.println(" \n");

            //  }
            long t2Final = System.nanoTime();
            System.out.println("time: "+(t2Final-t1Main)/1e9);
        }
    }
    }
    public static void doLJ(List<List<AtomType>> pairsAtoms, PotentialMoleculePair pTarget, LJUFF[] p2LJ, IPotential2[] p2lj, double rc, boolean doTrunc){
        PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
        int i = 0;
        UFF uff = new UFF();
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        //System.out.println(pairsAtoms + " Pairs");
        for(List<AtomType>individualPair: pairsAtoms){
            AtomType atomNameOne = individualPair.get(0);
            AtomType atomNameTwo = individualPair.get(1);
            String atomTypeStringOne = String.valueOf(atomNameOne);
            String atomTypeStringTwo = String.valueOf(atomNameTwo);
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);
            double[] iKey = pdbReaderReplica.atomicPot(atomTypeOne);
            double[] jKey = pdbReaderReplica.atomicPot(atomTypeTwo);
            double epsilonIKey = kcals.toSim(iKey[3]);
            double epsilonJKey = kcals.toSim(jKey[3]);
            double sigmaIKey = iKey[2];
            double sigmaJKey = jKey[2];
            //sigmaIJ[i] = (sigmaIKey + sigmaJKey) / 2;
            p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);
            if(doTrunc){
                TruncationFactory tf = new TruncationFactoryForceShift(rc);
                p2lj[i] = tf.make(p2LJ[i]);
                //System.out.println(atomTypeOne + " " + atomTypeTwo +" "+ Arrays.toString(iKey) +" " +Arrays.toString(jKey));
                pTarget.setAtomPotentials(atomNameOne, atomNameTwo, p2lj[i], new double[]{1, 0, 0, 1});
            }else {
                pTarget.setAtomPotentials(atomNameOne, atomNameTwo, p2LJ[i], new double[]{1, 0, 0, 1});
            }

            i++;
        }
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
                if(atomName1.equals("C_3") && atomName2.equals("O_1")){
                    bondOrder = 2;
                }
                if(atomName1.equals("O_3") && atomName2.equals("Cu")){
                    bondOrder = 1;
                }
                if(atomName1.equals("C_3") && atomName2.equals("H")){
                    bondOrder = 1;
                }
                if(atomName1.equals("C_3") && atomName2.equals("O_3")){
                    bondOrder = 1;
                }
                //  System.out.println(bondOrder +" bondorder");
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
                if(bondList1.size()< atom1){
                    bondList1.set(atom1, 1);
                } else if (bondList1.size()< atom2) {
                    bondList1.set(atom2, 1);
                }else if (bondList1.size()< atom3) {
                    bondList1.set(atom3, 1);
                }
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
        String confListName1 = "F://Avagadro//molecule//ethene";
        String confListName2 = "F://Avagadro//molecule//ethane";
        String confListName3 = "F://Avagadro//molecule//ch4";
        String confListName4 = "F://Avagadro//molecule//h2";
        String confListName5 = "F://Avagadro//molecule//Ar";

        String confName1 = "F://Avagadro//MOP_new//icosahedron//propyl_Co_MOP";
       // String confName2 = "F://Avagadro//mop//tetra_cu";
        String confName3 = "F://Avagadro//mop//tetra_cu_3C";
        String confName4 = "F://Avagadro//MOP_new//rhombic//propyl_Co_MOP";
      //  String confName2 = "F://Avagadro//mop//tetra_cu";
        // don't change these
        public int[] nTypes = new int[]{1, 1};
        public double temperature = 330;
        public  double truncLimit = 10.5;
        public double truncDiff = 5;
        public double start = 10;
        public long numSteps =100000;
        public double tempStart = 270;
        public double tempDiff = 30;
        public double tempFinal = 365;
        public Vector centreMoleculeOne = new Vector3D(0.0,0.0,0.0);
        public double refFrac = -1;
        public double sigmaHSRef = 9;
        //public Level level = Level.CLASSICAL;
        public boolean doHist = false;
        public boolean doTrunc = true;
      //  public Nonadditive nonAdditive = Nonadditive.NONE;
        public boolean useSZ = false;

    }

}
