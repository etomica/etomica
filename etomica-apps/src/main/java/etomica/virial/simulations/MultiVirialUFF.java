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
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential2;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.PotentialMoleculePair;
import etomica.potential.UFF.*;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.DimensionRatio;
import etomica.units.dimensions.Quantity;
import etomica.units.dimensions.Volume;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.MCMoveClusterAngleBend;
import etomica.virial.mcmove.MCMoveClusterStretch;
import etomica.virial.wheatley.ClusterWheatleySoftMix;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;

public class MultiVirialUFF {
    public static ISpecies species, species2;
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
        long t1Final = System.nanoTime();
        VirialUniversalParam params = new VirialUniversalParam();
        boolean isCommandline = args.length > 0;
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {

        }
        final int nPoints = params.nTypes[0] + params.nTypes[1];
        //   double temperature = params.temperature;
        long steps = params.numSteps;
        final int[] nTypes = params.nTypes;
        double refFreq = params.refFreq;
        //  double rc = params.rc;
        double sigmaHSRef = params.sigmaHSRef;
        Space space = Space3D.getInstance();
        List<String> confString = new ArrayList<>();
        confString.add(params.confName);
         confString.add(params.conf2Name);
        //  confString.add(params.conf3Name);
      //  for(int alpha =0; alpha< confString.size(); alpha++){
            String confName = confString.get(0);
            System.out.println(confName);
        String conf2Name = confString.get(1);
        System.out.println(conf2Name);
            PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
            species = pdbReaderReplica.getSpecies(confName, false,new Vector3D(10,0,0), false);
            //  System.out.println("Species");
            // box.addNewMolecule(species);
            PDBReaderMOP2 pdbReaderMOP2 = new PDBReaderMOP2();
            species2 = pdbReaderMOP2.getSpeciesMOP(conf2Name, false,new Vector3D(0,0,0), false);
            List<AtomType> atomTypes1 = species.getUniqueAtomTypes();
            List<AtomType> atomTypes2 = species2.getUniqueAtomTypes();
            List<List<AtomType>> pairsAtoms = new ArrayList<>();

            for (i = 0; i < atomTypes1.size(); i++){
                for (int j = 0; j < atomTypes2.size(); j++) {
                    if(i<=j){
                        List<AtomType> subPair = new ArrayList<>();
                        subPair.add(species.getAtomType(i));
                     //   System.out.println(species.getAtomType(i) + " " +species2.getAtomType(j)  + " Pair");
                        subPair.add(species2.getAtomType(j));
                        pairsAtoms.add(subPair);
                    //    System.out.println(subPair + " subPair");
                    }
                }
            }
            //System.exit(1);
            int atom1 , atom2 , atom3;
            String atomName1 , atomName2, atomName3;

            int pairAtomSize = pairsAtoms.size();
            // System.out.println(pairsAtoms + " pairs " + pairAtomSize);
            SpeciesBuilder speciesBuilder = new SpeciesBuilder(Space3D.getInstance());
            speciesBuilder.setDynamic(true);
            SpeciesManager sm1 = new SpeciesManager.Builder().addSpecies(species2).addSpecies(species).build();

            PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm1);
            // System.out.println(nPoints+" at "+temperature+"K");
            // temperature = Kelvin.UNIT.toSim(temperature);
            ArrayList<ArrayList<Integer>> connectedAtoms = pdbReaderReplica.getConnectivityWithoutRunning();
            ArrayList<ArrayList<Integer>> connectivityModified = pdbReaderReplica.getConnectivityModifiedWithoutRunning();
            Map<Integer,String> atomMap = pdbReaderReplica.getAtomMapWithoutRunning();
            HashMap<Integer, String> atomMapModified = pdbReaderReplica.getAtomMapModifiedWithoutRunning();
            List<int[]> duplets = pdbReaderReplica.getBondedAtomList(connectivityModified);
            List<int[]> triplets = pdbReaderReplica.getAngleList(connectivityModified);
            List<int[]> quadruplets = pdbReaderReplica.getTorsionList(connectedAtoms);
            ArrayList<Integer> bondList = pdbReaderReplica.getBondList(connectedAtoms, atomMap);
            //  System.out.println(bondList +" bondList");
            Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT),Mole.UNIT);
            Map<String, double[]> atomicPotMap = pdbReaderReplica.atomicPotMap();
            makeAtomPotentials(sm1);
            Map<Integer, String> atomIdentifierMapModified = pdbReaderReplica.atomIdentifierMapModified(connectivityModified, atomMapModified);
            //  System.out.println(atomIdentifierMapModified + "atomIdentifierMapModified");
            List<int[]>dupletsSorted= pdbReaderReplica.getDupletesSorted();
            List<int[]>tripletsSorted= pdbReaderReplica.getAnglesSorted();
            List<int[]>quadrupletsSorted= pdbReaderReplica.getTorsionSorted();
            //ArrayList<Integer> bondsNum = PDBReader.getBonds();
            Map<String[],List<int[]>> bondTypesMap= pdbReaderReplica.idenBondTypes(dupletsSorted, atomIdentifierMapModified);
            Map<String[],List<int[]>> angleTypesMap= pdbReaderReplica.idenAngleTypes(tripletsSorted, atomIdentifierMapModified);
            Map<String[],List<int[]>> torsionTypesMap= pdbReaderReplica.idenTorsionTypes(quadrupletsSorted, atomIdentifierMapModified);
            //   System.out.println(connectivityModified);
            // PDBReaderReplica.centreOfMass(species, atomIdentifierMapModified);

            ArrayList<ArrayList<Integer>> modifiedOutput = new ArrayList<>();
            for (ArrayList<Integer> innerList : connectivityModified) {
                ArrayList<Integer> modifiedInnerList = new ArrayList<>(innerList.subList(1, innerList.size()));
                modifiedOutput.add(modifiedInnerList);
            }
            //   System.out.println(modifiedOutput +" modified");
            IntArrayList[] dupletsIntArrayList = new IntArrayList[modifiedOutput.size()];

            for (i = 0; i < modifiedOutput.size(); i++) {
                ArrayList<Integer> innerList = modifiedOutput.get(i);
                IntArrayList intArrayList = new IntArrayList(innerList.size());
                for (int j = 0; j < innerList.size(); j++) {
                    intArrayList.add(innerList.get(j));
                    // System.out.println(intArrayList);
                }
                dupletsIntArrayList[i] = intArrayList;
                // System.out.println(dupletsIntArrayList[i]);
            }
      /*  for (IntArrayList list : dupletsIntArrayList) {
            for (i = 0; i < list.size(); i++) {
                int value = list.getInt(i);
                System.out.print(value + " ");
            }
            System.out.println();
        }*/

            Set<String> uniqueAtoms = pdbReaderReplica.uniqueElementIdentifier();
            // System.out.println(uniqueAtoms + "uniqueAtoms");
            //System.exit(0);
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

            PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm1);
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
                    //    System.out.println(Arrays.toString(bondParamsArray) + " ArrayToString");
                    bondConstant = UFF.BondConstantArray(bondParamsArray[0], bondParamsArray[1]);
                    //   System.out.println(Arrays.toString(bondConstant) + " arrayConstant");
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
                    int bondListValueTwo = 0;
                    if(atomName2.equals("C_1p") && atomName3.equals("O_1p")){
                        bondListValueTwo = 2;
                    }
                    int bondListValueOne = bondList.get(atom1);

                    int bondListValueThree = bondList.get(atom3);
                    double[] atomOnePot = atomicPotMap.get(atomName1);
                    double[] atomTwoPot = atomicPotMap.get(atomName2);
                    double[] atomThreePot = atomicPotMap.get(atomName3);
                    int num =0;
                    int caseNum =1;
                    angleParamsArray= UFF.angleUFF (atomOnePot[0], atomTwoPot[0], atomThreePot[0], atomOnePot[5], atomTwoPot[5], atomThreePot[5], atomOnePot[6], atomTwoPot[6],atomThreePot[6], atomTwoPot[1], bondListValueOne, bondListValueTwo, bondListValueThree,0);
                    //   System.out.println(Arrays.toString(angleParamsArray) + " arrayAngle");
                    P3BondAngleUFF p3Angle = new P3BondAngleUFF(angleParamsArray[0],  angleParamsArray[1], angleParamsArray[2], angleParamsArray[3],atomTwoPot[1], 0, caseNum);
                    bondingInfo.setBondingPotentialTriplet(species, p3Angle, angle);
                    break;
                }

            }

            //  System.out.println(species+"Species");
            P4BondTorsionUFF torsionUFF = null;
            double Vi =0, Vj =0, V=0, Vtrue=0,  type;
            int p;
            //    System.out.println(quadruplets.size() + " Quad size");

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
            //System.out.println(Arrays.deepToString(p4ValueArray.toArray()));
            // set pairwise potentials
            LJUFF[] p2LJ = new LJUFF[pairAtomSize];
            IPotential2[] p2lj = new IPotential2[pairAtomSize];
            double[] sigmaIJ = new double[pairAtomSize];
            i = 0;
            //  System.out.println(Arrays.deepToString(pairsAtoms.toArray()));
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
                double[] iKey = pdbReaderReplica.atomicPot(atomTypeOne);
                double[] jKey = pdbReaderReplica.atomicPot(atomTypeTwo);
                //System.out.println(atomNameOne + " " + Arrays.toString(iKey) + " " + atomNameTwo + " " + Arrays.toString(jKey));
                double epsilonIKey = Kelvin.UNIT.toSim(iKey[1]);
               // System.out.println(atomNameOne +" "+ atomNameTwo);
                double epsilonJKey = kcals.toSim(jKey[3]);
                double sigmaIKey = 1.12246 *  iKey[0];
                double sigmaJKey = jKey[2];
                sigmaIJ[i] = (sigmaIKey + sigmaJKey) / 2;
                //    TruncationFactoryForceShift tf = new TruncationFactoryForceShift(rc);
                p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);
                //  p2lj[i] = tf.make(p2LJ[i]);
                // p2lj[i] = uff.vdwNew(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey, tf);
              //  System.out.println(atomNameOne + " " + atomNameTwo);
                // potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo,p2lj[i]);
                pTargetAB.setAtomPotentials(atomNameOne, atomNameTwo, p2LJ[i], new double[]{1, 0, 0, 1});
                i++;
            }

        List<List<AtomType>> pairsAtoms2 = new ArrayList<>();

        for (i = 0; i < atomTypes1.size(); i++){
            for (int j = 0; j < atomTypes1.size(); j++) {
                if(i<=j){
                    List<AtomType> subPair = new ArrayList<>();
                    subPair.add(species.getAtomType(i));
                    //   System.out.println(species.getAtomType(i) + " " +species2.getAtomType(j)  + " Pair");
                    subPair.add(species.getAtomType(j));
                    pairsAtoms2.add(subPair);
                    //    System.out.println(subPair + " subPair");
                }
            }
        }
        for(List<AtomType> atomPairs : pairsAtoms2) {
            //String nameKey ="key";
            AtomType atomNameOne = atomPairs.get(0);
            AtomType atomNameTwo = atomPairs.get(1);
            // System.out.println( atomNameOne + " " + atomNameTwo);
            String atomTypeStringOne = String.valueOf(atomPairs.get(0));
            String atomTypeStringTwo = String.valueOf(atomPairs.get(1));
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);
            double[] iKey = pdbReaderReplica.atomicPot(atomTypeOne);
            double[] jKey = pdbReaderReplica.atomicPot(atomTypeTwo);
            //System.out.println(atomNameOne + " " + Arrays.toString(iKey) + " " + atomNameTwo + " " + Arrays.toString(jKey));
            //System.out.println(atomNameOne +" "+ atomNameTwo);
            double epsilonIKey = Kelvin.UNIT.toSim(iKey[1]);
            double epsilonJKey = Kelvin.UNIT.toSim(jKey[1]);
            double sigmaIKey = 1.12246 *  iKey[0];
            double sigmaJKey = 1.12246 *  jKey[0];
            sigmaIJ[i] = (sigmaIKey + sigmaJKey) / 2;
            //    TruncationFactoryForceShift tf = new TruncationFactoryForceShift(rc);
            p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);
            //  p2lj[i] = tf.make(p2LJ[i]);
            // p2lj[i] = uff.vdwNew(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey, tf);
           // System.out.println(atomNameOne + " " + atomNameTwo);
            // potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo,p2lj[i]);
            pTargetA.setAtomPotentials(atomNameOne, atomNameTwo, p2LJ[i], new double[]{1, 0, 0, 1});
            i++;
        }
       // System.exit(1);
            List<Integer> temp = new ArrayList<>();

           // for(int j =1; j< 2; j++){
                double temperatureNew =Kelvin.UNIT.toSim(params.temperature) ;
                System.out.println("\n" +params.temperature + " " + temperatureNew);
                targetCluster.setTemperature(temperatureNew);
                refCluster.setTemperature(temperatureNew);
                long newStep = params.numSteps;
                for (i=0; i<targetDiagrams.length; i++) {
                    targetDiagrams[i].setTemperature(temperatureNew);
                }

                System.out.println("sigmaHSRef: "+params.sigmaHSRef);
                // eovererr expects this string, BnHS
                System.out.println("B"+nPoints+"HS: "+refIntegral);
             //   final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, sm1,
                     //   new int[]{alkaneFlex ? (nPoints+1) : nPoints},temperatureNew, refCluster, targetCluster);
                final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species,species2}, nTypes, temperatureNew, refCluster, targetCluster);
                sim.setExtraTargetClusters(targetDiagrams);
                sim.setBondingInfo(bondingInfo);
                sim.setIntraPairPotentials(pTargetA.getAtomPotentials());
                // sim.setRandom(new RandomMersenneTwister(2));
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
                    //  ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
                    // ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[0]).setConstraintMap(constraintMap);
                    // ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).setConstraintMap(constraintMap);
                }

                if (refFreq >= 0) {
                    sim.integratorOS.setAdjustStepFraction(false);
                    sim.integratorOS.setRefStepFraction(refFreq);
                }

                sim.integratorOS.setNumSubSteps(1000);
                sim.integratorOS.setAggressiveAdjustStepFraction(true);
                System.out.println(newStep+" steps (100 blocks of "+newStep/1000+")");
                newStep /= 1000;



                sim.integratorOS.setNumSubSteps(1000);
                long t1 = System.nanoTime();
        if(false) {
            double size =20;
            sim.box[0].getBoundary().setBoxSize(etomica.space.Vector.of(new double[]{size, size, size}));
            sim.box[1].getBoundary().setBoxSize(etomica.space.Vector.of(new double[]{size, size, size}));
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
            //    IDataInfo dataInfo = new DataDouble.DataInfoDouble("B" + nPoints, new CompoundDimension(new Dimension[]{new Dimension(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints - 1}));
            //  Unit unit = new CompoundUnit(new Unit[]{new UnitRatio(Liter.UNIT, Mole.UNIT)}, new double[]{nPoints - 1});
            //  averageBox.putDataInfo(dataInfo);
            ////  averageBox.setLabel("average");
            // averageBox.setUnit(unit);
            //  errorBox.putDataInfo(dataInfo);
            //errorBox.setLabel("error");
            //errorBox.setPrecision(2);
            // errorBox.setUnit(unit);
            sim.integratorOS.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));
            return;
        }
                // if running interactively, don't use the file
                String refFileName = isCommandline ? "refpref"+nPoints+"_"+params.temperature : null;
                // this will either read the refpref in from a file or run a short simulation to find it
                sim.initRefPref(refFileName, steps/40);


                // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
                // if it does continue looking for a pref, it will write the value to the file
                sim.equilibrate(refFileName, steps/20);
                ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, 1000);
                System.out.println(steps);
                sim.setAccumulatorBlockSize(newStep);
                System.out.println("\n +" + sim.temperature + " K" );
                System.out.println("equilibration finished");
                sim.setAccumulatorBlockSize(newStep);
                sim.integratorOS.setNumSubSteps((int)newStep);
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
                System.out.println("\n");
       //     }

            long t2Final = System.nanoTime();
            System.out.println("time: "+(t2Final-t1Final)/1e9);





       // }
    }
    public static class VirialUniversalParam extends ParameterBase {
        public String confName ="Xe";
        public String conf2Name ="output_structure_min";
        public int[] nPoints = new int[]{1, 1};
        public double temperature =300;// Kelvin
        public long numSteps = 1000000;
        public double rc = 10;
        public double refFreq = -1;
        public double sigmaHSRef = 5;
        public int[] nTypes = new int[]{1, 1};

    }
    public static IPotential2[][] makeAtomPotentials(SpeciesManager sm) {
        // we could try to store the potentials more compactly, but it doesn't really matter
        ISpecies species = sm.getSpecies(sm.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        // System.out.println(lastTypeIndex + 1+ " "+lastTypeIndex + 1 + " lastTypeIndex" + " "+species.getAtomType(species.getUniqueAtomTypeCount() - 1));
        return new IPotential2[lastTypeIndex + 1][lastTypeIndex + 1];
    }
}
