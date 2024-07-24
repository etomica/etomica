package etomica.potential.UFF;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
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
import java.util.List;
import java.util.*;


public class UniversalSimulation extends Simulation {

    public PotentialCompute pcAgg;
    public IntegratorMC integratorMC;
    // public SpeciesGeneral species;
    public Box box;
    public MCMoveVolume mcMoveVolume;
    public MCMoveMolecule mcMoveMolecule;
    public ISpecies species;
    double molecularWeight=0;

    public UniversalSimulation(Space space, double density, int numMoleules, double temperature, String configFileName, double rc, double pressure) {

        super(space);
        int atom1 = 0, atom2 = 0, atom3=0;
       // String confName = "F://Avagadro//mopstrands//demo//3" ;
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        int num = 1;
        Vector centreMol = new Vector3D(0.0,0.0,0.0);
        species = pdbReaderMOP.getSpecies(configFileName, false, centreMol, false);
        System.out.println("Species");
        System.out.println(species.getMass());
        List<AtomType> atomTypes = species.getUniqueAtomTypes();
        List<List<AtomType>> pairsAtoms = new ArrayList<>();
        System.out.println(species);
        for(int i=0; i<atomTypes.size(); i++) {
            for (int j = 0; j < atomTypes.size(); j++) {
                if(i<=j){
                    List<AtomType> subPair = new ArrayList<>();
                    subPair.add(species.getAtomType(i));
                    subPair.add(species.getAtomType(j));
                    pairsAtoms.add(subPair);
                    System.out.println(subPair + " subPair");
                }
            }
        }
        final boolean isPureAtoms = false;
        BondingInfo bondingInfo = null;
        int pairAtomSize = pairsAtoms.size();
        System.out.println(pairsAtoms + " pairs " + pairsAtoms.size());
        SpeciesBuilder speciesBuilder = new SpeciesBuilder(Space3D.getInstance());
        speciesBuilder.setDynamic(true);
        String atomName1 = null, atomName2= null, atomName3= null;
        setRandom(new RandomMersenneTwister(1));
        addSpecies(species);
        box = this.makeBox();
        box.setNMolecules(species, numMoleules);
        new BoxInflate(box, space, density).actionPerformed();
        System.out.println(species);

        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);
        System.out.println("In simulation Universal");
        ArrayList<ArrayList<Integer>> connectedAtoms = pdbReaderMOP.getConnectivity();
        System.out.println(connectedAtoms+ ": connectedAtom");
        ArrayList<ArrayList<Integer>> connectivityModified = pdbReaderMOP.getConnectivityModifiedWithoutRunning();
        System.out.println(connectivityModified+ ": connectedAtomModified" );
        Map<Integer,String> atomMap = pdbReaderMOP.getAtomMapWithoutRunning();
        System.out.println(atomMap + ": atomMap");
        Map<Integer, String> atomMapModified = pdbReaderMOP.getAtomMapModifiedWithoutRunning();
        System.out.println(atomMapModified + ": atomMapModified");
        List<int[]> duplets = pdbReaderMOP.getduplets();
        System.out.println(Arrays.deepToString(duplets.toArray())+ ": listOfBonds");
        List<int[]> triplets = pdbReaderMOP.gettriplets();
        System.out.println(Arrays.deepToString(triplets.toArray())+ ": listOfAngleModified");
        List<int[]> quadruplets = pdbReaderMOP.getquadruplets();
        System.out.println(Arrays.deepToString(quadruplets.toArray())+ " listOfTorsionModified");
        ArrayList<Integer> bondList = pdbReaderMOP.getBondList(connectedAtoms, atomMap);
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        Map<Integer, String> atomIdentifierMapModified = pdbReaderMOP.getModifiedAtomIdentifierMap();
        Map<String, double[]> atomicPotMap = pdbReaderMOP.atomicPotMap();
        System.out.println(atomicPotMap +" atomic Pot");
        makeAtomPotentials(sm);
        System.out.println(atomicPotMap + "atomicPotMap");
        System.out.println(atomIdentifierMapModified);
       // System.exit(1);
        ArrayList<Integer> bondsNum = pdbReaderMOP.getBonds();
        System.out.println(bondsNum + " bondNums");
        List<int[]> dupletsSorted= pdbReaderMOP.getDupletesSorted();
        Map<String, Integer> anglemap = pdbReaderMOP.getangleMap(dupletsSorted, bondsNum);
        System.out.println(Arrays.deepToString(dupletsSorted.toArray()) + " duplets Sorted");
        System.out.println(Arrays.deepToString(duplets.toArray())+ ": listOfBonds");
        List<int[]>tripletsSorted= pdbReaderMOP.getAnglesSorted();
        System.out.println(Arrays.deepToString(tripletsSorted.toArray()) + "triplets");
        List<int[]>quadrupletsSorted= pdbReaderMOP.getTorsionSorted();
        List<int[]> inversionsSorted = pdbReaderMOP.getListofInversions();
        System.out.println(Arrays.deepToString(inversionsSorted.toArray()));
        System.out.println("Here 137");
        //System.exit(1);
      // List<int[]>inversionSorted=PDBReader.getListofInversions();
       // System.out.println(Arrays.deepToString(inversionSorted.toArray()));

        //List<int[]>valuePairs = PDBReader.firstLocationOfAtomPrinter(atomIdentifierMapModified);
        //System.out.println(Arrays.deepToString(valuePairs.toArray()) + "valuePairs");
        Map<String[],List<int[]>> angleTypesMap= pdbReaderMOP.idenAngleTypes(tripletsSorted, atomIdentifierMapModified);
        Map<String[],List<int[]>> torsionTypesMap= pdbReaderMOP.idenTorsionTypes(quadrupletsSorted, atomIdentifierMapModified);
        Map<Integer, Integer> angleTypeSorter = pdbReaderMOP.coordinationNumberDeterminer(connectivityModified, atomIdentifierMapModified);
       // System.out.println(angleTypeSorter);
        //Set<String> uniqueAtoms = PDBReader.uniqueElementIdentifier();
        //System.out.println(uniqueAtoms + "uniqueAtoms");
        //Map<String, AtomType> typeMapNew = PDBReader.getTypeMapNew();
       // molecularWeight = PDBReader.getMolecularWeight();
        //System.exit(1);
        int i =0 ;
        UFF uff = new UFF();
        List<double[]> bondParams = new ArrayList<>();
        for(i=0; i<dupletsSorted.size(); i++){
            List<int[]> bonds = new ArrayList<>();
            int[] bondIndividual = dupletsSorted.get(i);
            bonds.add(bondIndividual);
            double[] bondParamsArray = new double[2];
            double[] bondConstant = new double[2];
            atom1 = bondIndividual[0];
            atom2 = bondIndividual[1];
            atomName1 = atomIdentifierMapModified.get(atom1);
            atomName2 = atomIdentifierMapModified.get(atom2);
            double[] atomOnePot = atomicPotMap.get(atomName1);
            double[] atomTwoPot = atomicPotMap.get(atomName2);
            double bondOrder = bondsNum.get(i);
            if(atomName1.equals("C_Ar") && atomName2.equals("C_Ar")){
                bondOrder = 1.5;
            } else {
                bondOrder = 1;
            }
            System.out.println(bondOrder +" bondorder");
          //  System.out.println(Arrays.toString(dupletsSorted.get(i)) + " " + bondOrder+ " " + atomName1 + " " + atomName2+" "+ Arrays.toString(atomOnePot) +" " + Arrays.toString(atomTwoPot));
            bondParamsArray= UFF.bondUFF (atomOnePot[0],  atomTwoPot[0], atomOnePot[5],  atomTwoPot[5], atomOnePot[6], atomTwoPot[6], bondOrder);
            //System.out.println(Arrays.toString(bondParamsArray) + " ArrayToString");
            bondConstant = UFF.BondConstantArray(bondParamsArray[0], bondParamsArray[1]);
            System.out.println(Arrays.toString(bondConstant) + " arrayConstant");
            P2HarmonicUFF p2Bond = new P2HarmonicUFF(bondParamsArray[0],  bondParamsArray[1]);
            bondParams.add(bondConstant);
            pmBonding.setBondingPotentialPair(species, p2Bond, bonds );
        }
        for( i=0; i<anglemap.size(); i++){
            //System.out.println(Arrays.deepToString(anglemap.keySet().toArray()) +  " AngleMap");
        }
     //   System.exit(1);
        // Print the map using toString()
        /*for (String key : anglemap.keySet()) {
            System.out.print("Key: [");
            for (String numNew : key) {
                System.out.print(numNew + " ");
            }
            System.out.println("] Value: " + anglemap.get(key));
        }*/
       // System.exit(1);
        //System.out.println(Arrays.deepToString(bondParams.toArray()));
        System.out.println("Start of angle");
        int[] arrOne = new int[2];
        int[] arrTwo = new int[2];
        bondingInfo = pmBonding.getBondingInfo();
       /* for(i=0; i<tripletsSorted.size(); i++){
            System.out.println("\n");
            System.out.println(bondsNum);

            if(tripletsSorted.get(i)[1] >tripletsSorted.get(i)[0]  ){
                arrOne = new int[]{tripletsSorted.get(i)[0], tripletsSorted.get(i)[1]};
            } else {
                arrOne = new int[]{tripletsSorted.get(i)[1], tripletsSorted.get(i)[0]};
            }
            if(tripletsSorted.get(i)[2] >tripletsSorted.get(i)[1]  ){
                arrTwo = new int[]{tripletsSorted.get(i)[1], tripletsSorted.get(i)[2]};
            } else {
               arrTwo = new int[]{tripletsSorted.get(i)[2], tripletsSorted.get(i)[1]};
            }
            System.out.println(Arrays.toString(arrOne) + " " + Arrays.toString(arrTwo));
            int[] arrOne = new int[]{tripletsSorted.get(i)[1], tripletsSorted.get(i)[0]};
            int[] arrTwo = new int[]{tripletsSorted.get(i)[1], tripletsSorted.get(i)[2]};
            int boOne = anglemap.get(arrOne);
            int boTwo = anglemap.get(arrTwo);
            System.out.println(Arrays.deepToString(tripletsSorted.toArray()) + " " + Arrays.toString(arrOne) + "  " + Arrays.toString(arrTwo) + " " + boOne + " " + boTwo);
            System.out.println(Arrays.toString(tripletsSorted.get(i)));
        }*/
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
                double bondListValueOne = bondList.get(atom1);
                double bondListValueTwo = bondList.get(atom2);
                if(atomName1.equals("C_2") && atomName2.equals("C_3")|| atomName1.equals("C_3") && atomName2.equals("C_2")){
                    bondListValueOne = 1;
                }
                bondListValueTwo = bondList.get(atom2);
                if(atomName3.equals("C_2") && atomName2.equals("C_3") || atomName3.equals("C_3") && atomName2.equals("C_2")){
                    bondListValueTwo = 1;
                }
                double bondListValueThree = bondList.get(atom3);
                double[] atomOnePot = atomicPotMap.get(atomName1);
                double[] atomTwoPot = atomicPotMap.get(atomName2);
                double[] atomThreePot = atomicPotMap.get(atomName3);
                int num2 = angleTypeSorter.get(atom2);
                System.out.println(atomName1 + " " + atomName2 + " " + atomName3);
                int caseNum =0;
               // if((atomName1.equals("C_2") ||  atomName2.equals("C_2"))|| atomName3.equals("C_2")){
                if(atomName2.equals("C_2")){
                    caseNum = 1;
                }
                if(num2 != 0 ){
                    num = num2;
                   // System.out.println( num +" " + Arrays.toString(angleIndividual));
                } else {
                    num = 0;
                   // System.out.println( num +" " + Arrays.toString(angleIndividual));
                }

                if(atomName1.equals("C_Ar") && atomName2.equals("C_Ar")){
                    bondListValueOne = 1.5;
                    bondListValueTwo = 1.5;
                    num=1;
                } else if (atomName1.equals("H") && atomName2.equals("C_Ar")) {
                    bondListValueOne = 1;
                    bondListValueTwo = 1.5;
                    num=1;
                } else if (atomName1.equals("C_3") && atomName2.equals("C_Ar") ||atomName2.equals("C_3") && atomName1.equals("C_Ar") ) {
                    bondListValueOne = 1;
                    bondListValueTwo = 1.5;
                }


                if (atomName2.equals("C_Ar") && atomName3.equals("C_Ar")) {
                    bondListValueTwo = 1.5;
                    bondListValueThree = 1.5;
                    num=1;
                } else if(atomName2.equals("C_Ar") && atomName3.equals("H")  ){
                    bondListValueTwo = 1.5;
                    bondListValueThree = 1;
                    num=1;
                } else if (atomName1.equals("C_3") && atomName2.equals("C_Ar")) {
                    bondListValueTwo = 1.5;
                    bondListValueThree = 1;
                }
                System.out.println(bondListValueOne +" " +bondListValueTwo +" " +bondListValueThree);
                angleParamsArray= UFF.angleUFF (atomOnePot[0], atomTwoPot[0], atomThreePot[0], atomOnePot[5], atomTwoPot[5], atomThreePot[5], atomOnePot[6], atomTwoPot[6],atomThreePot[6], atomTwoPot[1], bondListValueOne, bondListValueTwo, bondListValueThree, num);
                System.out.println(Arrays.toString(angleParamsArray) + " arrayAngle");
                double angleTwoPot = Degree.UNIT.toSim(atomTwoPot[1]);
                P3BondAngleUFF p3Angle = new P3BondAngleUFF(angleParamsArray[0],  angleParamsArray[1], angleParamsArray[2], angleParamsArray[3],angleTwoPot  ,num, caseNum);
                pmBonding.setBondingPotentialTriplet(species, p3Angle, angle);
                break;
            }

        }
       //System.exit(1);

      //  System.out.println(species+"Species");
        P4BondTorsionUFF torsionUFF = null;
        double Vi =0, Vj =0, V=0, Vtrue=0,  type;
        int p;
        System.out.println(quadruplets.size() + " Quad size");
        boolean quadNum, invertNum;
        if (quadruplets.size() ==0){
            quadNum = false;
           //System.out.println(quadNum + " quad "  + quadruplets.size() + " quad");
        } else {
            quadNum = true;
           // System.out.println(quadNum + " quad "+ quadruplets.size() + " quad");
        }

       if(inversionsSorted.size() ==0){
            invertNum= false;
        }  else {
            invertNum = true;

        }

       // System.out.println(invertNum + " quad "+ inversionSorted.size() + " quad");
        if(quadNum){
            for (Map.Entry<String[], List<int[]>> entry : torsionTypesMap.entrySet()) {
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
                    double bondOrder = 0;
                    if(atomName2.equals("C_Ar") && atomName3.equals("C_Ar")){
                         bondOrder = 1.5;
                    } else {
                         bondOrder = 1;
                    }
                    torsionParamsArray = UFF.torsionUFF(Vtrue, p, bondOrder);
                    P4BondTorsionUFF p4Torsion = new P4BondTorsionUFF(torsionParamsArray[0],(int) torsionParamsArray[1], torsionParamsArray[2]);
                    pmBonding.setBondingPotentialQuad(species, p4Torsion, torsion);
                    break;
                }

            }
        }
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 5, pmBonding.getBondingInfo());

        if(invertNum){
            int[] individualInversion = inversionsSorted.get(0);
            String atomOne = atomIdentifierMapModified.get(individualInversion[0]);
            if(atomOne.equals("C_2") || atomOne.equals("C_R")){
                int c0 = 1;
                int c1 =-1;
                int c2 = 0;
                double kijkl = kcals.toSim(6.0);
                P4BondInversionUFF p4Inversion = new P4BondInversionUFF(c0, c1, c2, kijkl);
                pmBonding.setBondingPotentialInvert(species, p4Inversion, inversionsSorted);
            } else{
                int c0 = 1;
                int c1 =-1;
                int c2 = 0;
                double kijkl = kcals.toSim(6.0);
                P4BondInversionUFF p4Inversion = new P4BondInversionUFF(c0, c1, c2, kijkl);
                pmBonding.setBondingPotentialInvert(species, p4Inversion, inversionsSorted);
            }
        }

        potentialMasterCell.doAllTruncationCorrection = false;
        LJUFF[] p2LJ = new LJUFF[pairAtomSize];
        IPotential2[] p2lj = new IPotential2[pairAtomSize];
        double[] sigmaIJ = new double[pairAtomSize];
       // System.out.println(pairsAtoms + " pairatoms");
        i = 0;
        System.out.println(Arrays.deepToString(pairsAtoms.toArray()));

        for(List<AtomType> atomPairs : pairsAtoms) {
            //String nameKey ="key";
            AtomType atomNameOne = atomPairs.get(0);
            AtomType atomNameTwo = atomPairs.get(1);
            // System.out.println( atomNameOne + " " + atomNameTwo);
            String atomTypeStringOne = String.valueOf(atomPairs.get(0));
            String atomTypeStringTwo = String.valueOf(atomPairs.get(1));
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);
            double[] iKey = pdbReaderMOP.atomicPot(atomTypeOne);
            double[] jKey = pdbReaderMOP.atomicPot(atomTypeTwo);
            //System.out.println(atomNameOne + " " + Arrays.toString(iKey) + " " + atomNameTwo + " " + Arrays.toString(jKey));
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


       // if (configFileName != null) {
           // ConfigurationFile config = new ConfigurationFile(configFileName);
           // config.initializeCoordinates(box);
           //BoxImposePbc.imposePBC(box);
        //}
      //  else {
            ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
            configuration.initializeCoordinates(box);
            potentialMasterCell.init();
            double u0 = potentialMasterCell.computeAll(false);
            double x = 1;
            System.out.println( u0 + " "+ x +" inMain "  + kcals.fromSim(u0));
            System.exit(1);
           while (u0 > 1e6*numMoleules) {
               //System.out.println(x +" before");
               x *= 0.99;
              // System.out.println( x +" =x");
                for(int j = 0; j< pairAtomSize; j++){
                    System.out.println(x + " " + sigmaIJ[0]);
                    p2LJ[j].setSigma(x*sigmaIJ[j]);
                    ((P2SoftSphericalSumTruncatedForceShifted)p2lj[j]).setTruncationRadius(rc);
                }
                u0 = potentialMasterCell.computeAll(false);
                System.out.println( u0 + " inMain afterwards " +kcals.fromSim(u0));
            }
            //System.out.println( u0 + " inMain afterwards "  + kcals.fromSim(u0));

            integratorMC.reset();
            while (u0 > 1e4*numMoleules) {
                while (u0 > 1e4 * numMoleules) {

                    integratorMC.doStep();
                    u0 = integratorMC.getPotentialEnergy();
                    // System.out.println("Inside Loop Two - 1");
                }
                while (x < 1 && u0 <= 1e4 * numMoleules) {
                    //  System.out.println("Inside Loop Two - 2");
                    x /= 0.99;
                    if (x > 1) x = 1;
                    for(int j = 0; j< pairAtomSize; j++){
                        //  System.out.println("Inside Loop Two - 3");
                        p2LJ[j].setSigma(x*sigmaIJ[j]);
                        ((P2SoftSphericalSumTruncatedForceShifted)p2lj[j]).setTruncationRadius(rc);
                    }
                    u0 = potentialMasterCell.computeAll(false);
                    //System.out.println(u0 +" inside Array @");
                }
                integratorMC.reset();
           // }
        }
    }

    public static void main(String[] args) throws IOException {
        final String APP_NAME = "Methane Universal";
        UniversalParams params = new UniversalParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.numSteps = 3500000;
            params.density = 0.00003;
            //params.configFilename = null; // "octane";
            params.graphics = true;
        }


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

        final UniversalSimulation sim = new UniversalSimulation(Space3D.getInstance(), density, numMolecules, temperature, configFilename, rc, pressure);

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

        if (false) {
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
            dhbt.setDiameter(sim.species.getTypeByName("C_2"), 3.54);
            //dhbt.setDiameter(sim.species.getTypeByName("C_3"), 3.54);
            dhbt.setDiameter(sim.species.getTypeByName("H"), 2.25);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getTypeByName("C_2"), Color.WHITE);
            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getTypeByName("C_3"), Color.WHITE);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getTypeByName("H"), Color.RED);

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
        try{
            writer.write("U: " +" "+avgU+ "   err: "+" "+errU+"   cor: " +" "+corU+ "\n");
            writer.write("P: " +" "+avgP+ "   err: " +" "+errP+ "   cor: " +corP+ "\n");
            writer.write("(Z-1)/rho: "+" "+avgZ_+"   err: "+" "+errZ_+"   cor: "+" "+corZ_+ "\n");
        } catch (IOException e) {
            System.out.println("An error occurred while writing to the file: " + e.getMessage());
        }
        writer.close();
        String absolutePath = file.getAbsolutePath();
        System.out.println("File path: " + absolutePath);
        System.out.println("time: "+" "+(t2-t1)/1e9);
    }

    public static class UniversalParams extends ParameterBase {
        public double temperatureK = 300;
        public int numMolecules = 1;
        //public int pressure = 10;
        public double density = 0.0000005;
        public boolean graphics = false;
        public long numSteps = 1000000;
        public String configFilename = "F://Avagadro//molecule//co2";
        public double rc = 10;
        public double pressureKPa = 1402;
    }


    public static IPotential2[][] makeAtomPotentials(SpeciesManager sm) {
        // we could try to store the potentials more compactly, but it doesn't really matter
        ISpecies species = sm.getSpecies(sm.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
       // System.out.println(lastTypeIndex + 1+ " "+lastTypeIndex + 1 + " lastTypeIndex" + " "+species.getAtomType(species.getUniqueAtomTypeCount() - 1));
        return new IPotential2[lastTypeIndex + 1][lastTypeIndex + 1];
    }
}

