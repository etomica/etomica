package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.graph.model.Graph;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.potential.UFF.*;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.*;
import etomica.units.*;
import etomica.units.dimensions.*;
import etomica.units.dimensions.Dimension;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.cluster.ClusterChainHS;
import etomica.virial.cluster.ClusterSum;
import etomica.virial.cluster.ClusterSumShell;
import etomica.virial.cluster.VirialDiagrams;
import etomica.virial.mcmove.*;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;

public class virialUFF {
    public static ISpecies species;
    public static Box box;

    public static String getSplitGraphString(Set<Graph> gSplit, VirialDiagrams flexDiagrams, boolean correction) {
        DeleteEdge edgeDeleter = new DeleteEdge();
        DeleteEdgeParameters ede = new DeleteEdgeParameters(flexDiagrams.eBond);
        boolean first = true;
        String str = "";
        for (Graph gs : gSplit) {
            if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.efbcBond)) {
                str += " "+gs.nodeCount()+"bc";
            }
            else {
                str += " "+gs.getStore().toNumberString();
                if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.eBond)) {
                    str += "p" + edgeDeleter.apply(gs, ede).getStore().toNumberString();
                }
            }
            if (first && correction) str += "c";
            first = false;
        }
        return str;
    }
    public static void main(String[] args) {
        int i =0;
        VirialUniversalParam params = new VirialUniversalParam();
        boolean isCommandline = args.length > 0;
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {

        }
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double refFreq = params.refFreq;
        double rc = params.rc;
        double sigmaHSRef = params.sigmaHSRef;
        Space space = Space3D.getInstance();
        String confName = params.confName;
        species = PDBReaderReplica.getSpecies(confName);
        System.out.println("Species");
       // box.addNewMolecule(species);

        java.util.List<AtomType> atomTypes = species.getUniqueAtomTypes();
        java.util.List<java.util.List<AtomType>> pairsAtoms = new ArrayList<>();
        System.out.println(species);
        for(i=0; i<atomTypes.size(); i++) {
            for (int j = 0; j < atomTypes.size(); j++) {
                if(i<=j){
                    List<AtomType> subPair = new ArrayList<>();
                    subPair.add(species.getAtomType(i));
                    System.out.println(species.getAtomType(i) + " " +species.getAtomType(j)  + " Pair");
                    subPair.add(species.getAtomType(j));
                    pairsAtoms.add(subPair);
                    //System.out.println(subPair + " subPair");
                }
            }
        }
        int atom1 , atom2 , atom3;
        String atomName1 , atomName2, atomName3;

        int pairAtomSize = pairsAtoms.size();
        System.out.println(pairsAtoms + " pairs " + pairAtomSize);
        SpeciesBuilder speciesBuilder = new SpeciesBuilder(Space3D.getInstance());
        speciesBuilder.setDynamic(true);
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm);
        System.out.println(nPoints+" at "+temperature+"K");
        temperature = Kelvin.UNIT.toSim(temperature);
        ArrayList<ArrayList<Integer>> connectedAtoms = PDBReaderReplica.getConnectivityWithoutRunning();
        ArrayList<ArrayList<Integer>> connectivityModified = PDBReaderReplica.getConnectivityModifiedWithoutRunning();
        Map<Integer,String> atomMap = PDBReaderReplica.getAtomMapWithoutRunning();
        HashMap<Integer, String> atomMapModified = PDBReaderReplica.getAtomMapModifiedWithoutRunning();
        List<int[]> duplets = PDBReaderReplica.getBondedAtomList(connectivityModified);
        List<int[]> triplets = PDBReaderReplica.getAngleList(connectivityModified);
        List<int[]> quadruplets = PDBReaderReplica.getTorsionList(connectedAtoms);
        ArrayList<Integer> bondList = PDBReaderReplica.getBondList(connectedAtoms, atomMap);
        System.out.println(bondList +" bondList");
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        Map<String, double[]> atomicPotMap = PDBReaderReplica.atomicPotMap();
        UniversalSimulation.makeAtomPotentials(sm);
        Map<Integer, String> atomIdentifierMapModified = PDBReaderReplica.atomIdentifierMapModified(connectivityModified, atomMapModified);
        System.out.println(atomIdentifierMapModified + "atomIdentifierMapModified");
        List<int[]>dupletsSorted= PDBReaderReplica.getDupletesSorted();
        List<int[]>tripletsSorted= PDBReaderReplica.getAnglesSorted();
        List<int[]>quadrupletsSorted= PDBReaderReplica.getTorsionSorted();
        //ArrayList<Integer> bondsNum = PDBReader.getBonds();
        Map<String[],List<int[]>> bondTypesMap= PDBReaderReplica.idenBondTypes(dupletsSorted, atomIdentifierMapModified);
        Map<String[],List<int[]>> angleTypesMap= PDBReaderReplica.idenAngleTypes(tripletsSorted, atomIdentifierMapModified);
        Map<String[],List<int[]>> torsionTypesMap= PDBReaderReplica.idenTorsionTypes(quadrupletsSorted, atomIdentifierMapModified);
        System.out.println(connectivityModified);
       // PDBReaderReplica.centreOfMass(species, atomIdentifierMapModified);

        ArrayList<ArrayList<Integer>> modifiedOutput = new ArrayList<>();
        for (ArrayList<Integer> innerList : connectivityModified) {
            ArrayList<Integer> modifiedInnerList = new ArrayList<>(innerList.subList(1, innerList.size()));
            modifiedOutput.add(modifiedInnerList);
        }
        System.out.println(modifiedOutput +" modified");
        IntArrayList[] dupletsIntArrayList = new IntArrayList[modifiedOutput.size()];

        for (i = 0; i < modifiedOutput.size(); i++) {
            ArrayList<Integer> innerList = modifiedOutput.get(i);
            IntArrayList intArrayList = new IntArrayList(innerList.size());
            for (int j = 0; j < innerList.size(); j++) {
                intArrayList.add(innerList.get(j));
                System.out.println(intArrayList);
            }
            dupletsIntArrayList[i] = intArrayList;
            System.out.println(dupletsIntArrayList[i]);
        }
        for (IntArrayList list : dupletsIntArrayList) {
            for (i = 0; i < list.size(); i++) {
                int value = list.getInt(i);
                System.out.print(value + " ");
            }
            System.out.println();
        }

        Set<String> uniqueAtoms = PDBReaderReplica.uniqueElementIdentifier();
        System.out.println(uniqueAtoms + "uniqueAtoms");
        //System.exit(0);
        box = new Box(space);
        // FIll it up with Meyer Function and Target Diagram
        MayerGeneral fTarget = new MayerGeneral(pTarget);
        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }
            @Override
            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef * sigmaHSRef ? 1 : 0;
            }
        };

        boolean alkaneFlex = nPoints > 2;
        VirialDiagrams alkaneDiagrams = new VirialDiagrams(nPoints, false, alkaneFlex);
        alkaneDiagrams.setDoReeHoover(false);
        ClusterSum targetCluster = alkaneDiagrams.makeVirialCluster(fTarget);
        ClusterChainHS refCluster = new ClusterChainHS(nPoints, fRefPos);
        double vhs = (4.0 / 3.0) * Math.PI * params.sigmaHSRef * params.sigmaHSRef *params.sigmaHSRef;
        final double refIntegral = SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1);

        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
        int[] targetDiagramNumbers = new int[0];
        boolean[] diagramFlexCorrection = null;
        targetDiagrams = alkaneDiagrams.makeSingleVirialClusters(targetCluster, null, fTarget);
        targetDiagramNumbers = new int[targetDiagrams.length];
        System.out.println("individual clusters:");
        Set<Graph> singleGraphs = alkaneDiagrams.getMSMCGraphs(true, false);
        Map<Graph,Graph> cancelMap = alkaneDiagrams.getCancelMap();
        int iGraph = 0;
        diagramFlexCorrection = new boolean[targetDiagrams.length];
        for (Graph g : singleGraphs) {
            System.out.print(iGraph+" ("+g.coefficient()+") "+g.getStore().toNumberString()); // toNumberString: its corresponding number
            targetDiagramNumbers[iGraph] = Integer.parseInt(g.getStore().toNumberString());

            Graph cancelGraph = cancelMap.get(g);
            if (cancelGraph != null) {
                diagramFlexCorrection[iGraph] = true;
                Set<Graph> gSplit = alkaneDiagrams.getSplitDisconnectedVirialGraphs(cancelGraph);

                System.out.print(" - "+getSplitGraphString(gSplit, alkaneDiagrams, false));

            }
            System.out.println();
            iGraph++;
        }


        System.out.println();
        Set<Graph> disconnectedGraphs = alkaneDiagrams.getExtraDisconnectedVirialGraphs();
        if (disconnectedGraphs.size() > 0) {
            System.out.println("extra clusters:");

            for (Graph g : disconnectedGraphs) {
                Set<Graph> gSplit = alkaneDiagrams.getSplitDisconnectedVirialGraphs(g);
                System.out.println(g.coefficient()+" "+getSplitGraphString(gSplit, alkaneDiagrams, true));
            }
            System.out.println();
        }


        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        for (i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }

        System.out.println("sigmaHSRef: "+params.sigmaHSRef);
        // eovererr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm);
        for (Map.Entry<String[], List<int[]>> entry : bondTypesMap.entrySet()) {
            String[] bondType = entry.getKey();
            List<int[]> bonds = entry.getValue();
            // System.out.println(Arrays.toString(bondType) + ": " + Arrays.deepToString(bonds.toArray()));
            for(int[]bondIndividual: bonds){

                double[] bondParamsArray = new double[2];
                double[] bondConstant = new double[2];
                atom1 = bondIndividual[0];
                atom2 = bondIndividual[1];
                atomName1 = atomIdentifierMapModified.get(atom1);
                atomName2 = atomIdentifierMapModified.get(atom2);
                double[] atomOnePot = atomicPotMap.get(atomName1);
                double[] atomTwoPot = atomicPotMap.get(atomName2);
                int bondListValueOne = bondList.get(atom1);
                int bondListValueTwo = bondList.get(atom2);
                int bondOrder = 1;
                bondParamsArray= UFF.bondUFF(atomOnePot[0],  atomTwoPot[0], atomOnePot[5],  atomTwoPot[5], atomOnePot[6], atomTwoPot[6] , bondOrder);
                System.out.println(Arrays.toString(bondParamsArray) + " ArrayToString");
                bondConstant = UFF.BondConstantArray(bondParamsArray[0], bondParamsArray[1]);
                System.out.println(Arrays.toString(bondConstant) + " arrayConstant");
                P2HarmonicUFF p2Bond = new P2HarmonicUFF(bondParamsArray[0],  bondParamsArray[1]);
                bondingInfo.setBondingPotentialPair(species, p2Bond, bonds);
                i++;
                break;
            }

        }

        for (Map.Entry<String[], List<int[]>> entry : angleTypesMap.entrySet()) {
            String[] angleType = entry.getKey();
            List<int[]> angle = entry.getValue();
            // System.out.println(Arrays.toString(angleType) + ": " + Arrays.deepToString(angle.toArray()));
            for(int[]angleIndividual: angle){
                double[] angleParamsArray = new double[4];
                atom1 = angleIndividual[0];
                atom2 = angleIndividual[1];
                atom3 = angleIndividual[2];
                atomName1 = atomIdentifierMapModified.get(atom1);
                atomName2 = atomIdentifierMapModified.get(atom2);
                atomName3 = atomIdentifierMapModified.get(atom3);
                int bondListValueOne = bondList.get(atom1);
                int bondListValueTwo = bondList.get(atom2);
                int bondListValueThree = bondList.get(atom3);
                double[] atomOnePot = atomicPotMap.get(atomName1);
                double[] atomTwoPot = atomicPotMap.get(atomName2);
                double[] atomThreePot = atomicPotMap.get(atomName3);
                int num =0;
                int caseNum =1;
                angleParamsArray= UFF.angleUFF (atomOnePot[0], atomTwoPot[0], atomThreePot[0], atomOnePot[5], atomTwoPot[5], atomThreePot[5], atomOnePot[6], atomTwoPot[6],atomThreePot[6], atomTwoPot[1], bondListValueOne, bondListValueTwo, bondListValueThree,0);
                System.out.println(Arrays.toString(angleParamsArray) + " arrayAngle");
                P3BondAngleUFF p3Angle = new P3BondAngleUFF(angleParamsArray[0],  angleParamsArray[1], angleParamsArray[2], angleParamsArray[3],atomTwoPot[1], 0, caseNum);
                bondingInfo.setBondingPotentialTriplet(species, p3Angle, angle);
                break;
            }

        }

        System.out.println(species+"Species");
        P4BondTorsionUFF torsionUFF = null;
        double Vi =0, Vj =0, V=0, Vtrue=0,  type;
        int p;
        System.out.println(quadruplets.size() + " Quad size");

        P4BondTorsionUFF[] p4BondTorsionArray = new P4BondTorsionUFF[quadruplets.size()];
        ArrayList<double[]> p4ValueArray = new ArrayList<>();

        for (Map.Entry<String[], List<int[]>> entry : torsionTypesMap.entrySet()) {
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
                atomName2 = atomIdentifierMapModified.get(atom2);
                atomName3 = atomIdentifierMapModified.get(atom3);
                Vi = UFF.switchCaseTorsion(atomName2);
                int bondListValueOne = bondList.get(torsionIndividual[1]);
                int bondListValueTwo = bondList.get(torsionIndividual[2]);
                p = p + bondListValueOne + bondListValueTwo;
                Vj = UFF.switchCaseTorsion(atomName3);
                V = Math.sqrt(Vi*Vj);
                Vtrue = kcals.toSim(V);
                double bondOrder = 1 ;
                torsionParamsArray = UFF.torsionUFF(Vtrue, p, bondOrder);
                p4BondTorsionArray[i] = new P4BondTorsionUFF(torsionParamsArray[0], (int) torsionParamsArray[1], torsionParamsArray[2]);
                double[] array = {torsionParamsArray[0], torsionParamsArray[1], torsionParamsArray[2]};
                p4ValueArray.add(array);
                bondingInfo.setBondingPotentialQuad(species, p4BondTorsionArray[i], torsion);
                i++;
            }
        }
        System.out.println(Arrays.deepToString(p4ValueArray.toArray()));

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, sm,
                new int[]{alkaneFlex ? (nPoints+1) : nPoints},temperature, refCluster, targetCluster);
        sim.setExtraTargetClusters(targetDiagrams);
        sim.setBondingInfo(bondingInfo);
        sim.setIntraPairPotentials(pTarget.getAtomPotentials());
        sim.setRandom(new RandomMersenneTwister(2));
        sim.init();
        PotentialCompute pc0 = sim.integrators[0].getPotentialCompute();
        //MCMoveClusterMoleculeMulti mcMoveMulti = new MCMoveClusterMoleculeMulti(pc0, space, sim.getRandom(), 0.05);
        if(bondTypesMap.size()> 0){
            MCMoveClusterStretch mcMoveStretchRef = new MCMoveClusterStretch(pc0, space, dupletsIntArrayList, sim.getRandom(), 0.05);
            MCMoveClusterAngleBend mcMoveBendRef = new MCMoveClusterAngleBend(pc0, sim.getRandom(), 0.05, space);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveStretchRef);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveBendRef);
            PotentialCompute pc1 = sim.integrators[1].getPotentialCompute();
            MCMoveClusterStretch mcMoveStretchTarget = new MCMoveClusterStretch(pc1, space, dupletsIntArrayList, sim.getRandom(), 0.05);
            MCMoveClusterAngleBend mcMoveBendTarget = new MCMoveClusterAngleBend(pc1, sim.getRandom(), 0.05, space);
            sim.integrators[1].getMoveManager().addMCMove(mcMoveStretchTarget);
            sim.integrators[1].getMoveManager().addMCMove(mcMoveBendTarget);
        } //else {
           // MCMoveClusterMoleculeMulti mcMoveClusterMoleculeRef = new MCMoveClusterMoleculeMulti(sim.getRandom(), box);
      //  }

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
        System.out.println(steps+" steps (100 blocks of "+steps/1000+")");
        steps /= 1000;

        // set pairwise potentials
        LJUFF[] p2LJ = new LJUFF[pairAtomSize];
        IPotential2[] p2lj = new IPotential2[pairAtomSize];
        double[] sigmaIJ = new double[pairAtomSize];
        i = 0;
        System.out.println(Arrays.deepToString(pairsAtoms.toArray()));
        UFF uff = new UFF();
        for(List<AtomType> atomPairs : pairsAtoms) {
            //String nameKey ="key";
            AtomType atomNameOne = atomPairs.get(0);
            AtomType atomNameTwo = atomPairs.get(1);
           // System.out.println( atomNameOne + " " + atomNameTwo);
            String atomTypeStringOne = String.valueOf(atomPairs.get(0));
            String atomTypeStringTwo = String.valueOf(atomPairs.get(1));
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);
            double[] iKey = PDBReaderReplica.atomicPot(atomTypeOne);
            double[] jKey = PDBReaderReplica.atomicPot(atomTypeTwo);
            //System.out.println(atomNameOne + " " + Arrays.toString(iKey) + " " + atomNameTwo + " " + Arrays.toString(jKey));
            double epsilonIKey = kcals.toSim(iKey[3]);
            double epsilonJKey = kcals.toSim(jKey[3]);
            double sigmaIKey = iKey[2];
            double sigmaJKey = jKey[2];
            sigmaIJ[i] = (sigmaIKey + sigmaJKey) / 2;
            TruncationFactory tf = new TruncationFactoryForceShift(rc);
            p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);
            p2lj[i] = tf.make(p2LJ[i]);
            System.out.println(atomNameOne + " " + atomNameTwo);
            // potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo,p2lj[i]);
            pTarget.setAtomPotentials(atomNameOne, atomNameTwo, p2lj[i], new double[]{1, 0, 0, 1});
            i++;
        }

        sim.integratorOS.setNumSubSteps(1000);


        if(false) {
            double size =20;
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
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
    public static class VirialUniversalParam extends ParameterBase {
        public String confName ="F://Avagadro//molecule//butane";
        public int nPoints =4;
        public double temperature = 100;// Kelvin
        public long numSteps = 500000;
        public double refFreq = -1;
        public double sigmaHSRef = 5.23;
        public double rc = 10;
    }
}
