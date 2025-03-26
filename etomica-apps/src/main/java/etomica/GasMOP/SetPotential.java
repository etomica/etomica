package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.potential.GAFF.GAFFSimulation;
import etomica.potential.OPLS_AA.PDBDataExtracterOPLS;
import etomica.potential.TraPPE.SpeciesGasTraPPE;
import etomica.potential.UFF.*;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.*;

import java.util.*;

public class SetPotential {
    public IntegratorMC integrator;
  //  public MCMoveInsertDelete mcMoveID;
    public Box box;
    public P2LennardJones potential;
    public SpeciesManager sm;
    public static int atom1, atom2, atom3;
    public static double sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey;
    public int i;
    public static String atomName1,atomName2, atomName3 ;


    public void setBondStretch(ISpecies species1, Map<String[], List<int[]>> bondTypesMap1, Map<String[],List<int[]>> angleTypesMap1, Map<String[],List<int[]>> torsionTypesMap1, ArrayList<Integer> bondsNum1, ArrayList<Integer> bondList1, List<int[]>quadrupletsSorted1, Map<Integer, String> atomIdentifierMapModified1, Map<String, double[]> atomicPotMap1, PotentialMasterBonding pmBonding){
        double Vi =0, Vj =0, V=0, Vtrue=0,  type;
        int p;
        int i =0;
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT),Mole.UNIT);
        for (Map.Entry<String[], List<int[]>> entry : bondTypesMap1.entrySet()) {
            String[] bondType = entry.getKey();
            List<int[]> bonds = entry.getValue();
            // System.out.println(Arrays.toString(bondType) + ": " + Arrays.deepToString(bonds.toArray()));
            for(int[]bondIndividual: bonds){
              //  bonds.add(bondIndividual);
                double[] bondParamsArray = new double[2];
                double[] bondConstant = new double[2];
                atom1 = bondIndividual[0];
                atom2 = bondIndividual[1];
                atomName1 = atomIdentifierMapModified1.get(atom1);
                atomName2 = atomIdentifierMapModified1.get(atom2);
                double[] atomOnePot = atomicPotMap1.get(atomName1);
                double[] atomTwoPot = atomicPotMap1.get(atomName2);
                double bondOrder = bondsNum1.get(i);
                if(atomName1.equals("C_Ar") && atomName2.equals("C_Ar")){
                    bondOrder = 1.5;
                }
                if(atomName1.equals("C_1") && atomName2.equals("O_1")){
                    bondOrder =2;
                }
                //System.out.println(bondOrder +" bondorder");
             //     System.out.println(Arrays.toString(dupletsSorted.get(i)) + " " + bondOrder+ " " + atomName1 + " " + atomName2+" "+ Arrays.toString(atomOnePot) +" " + Arrays.toString(atomTwoPot));
                bondParamsArray= UFF.bondUFF (atomOnePot[0],  atomTwoPot[0], atomOnePot[5],  atomTwoPot[5], atomOnePot[6], atomTwoPot[6], bondOrder);
             //   System.out.println(Arrays.toString(bondParamsArray) + " ArrayToString");
                // bondConstant = UFF.BondConstantArray(bondParamsArray[0], bondParamsArray[1]);
                //System.out.println(Arrays.toString(bondConstant) + " arrayConstant");
                P2HarmonicUFF p2Bond = new P2HarmonicUFF(bondParamsArray[0],  bondParamsArray[1]);
                //bondParams.add(bondConstant);
                pmBonding.setBondingPotentialPair(species1, p2Bond, bonds );
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
                int bondListValueOne = 1;
                int bondListValueTwo = 1;
                int bondListValueThree = 1;
                double[] atomOnePot = atomicPotMap1.get(atomName1);
                double[] atomTwoPot = atomicPotMap1.get(atomName2);
                double[] atomThreePot = atomicPotMap1.get(atomName3);
                int num =0;
                int caseNum =1;
                angleParamsArray= UFF.angleUFF (atomOnePot[0], atomTwoPot[0], atomThreePot[0], atomOnePot[5], atomTwoPot[5], atomThreePot[5], atomOnePot[6], atomTwoPot[6],atomThreePot[6], atomTwoPot[1], bondListValueOne, bondListValueTwo, bondListValueThree,0);
                // System.out.println(Arrays.toString(angleParamsArray) + " arrayAngle");
                P3BondAngleUFF p3Angle = new P3BondAngleUFF(angleParamsArray[0],  angleParamsArray[1], angleParamsArray[2], angleParamsArray[3],atomTwoPot[1], 0, caseNum);
                pmBonding.setBondingPotentialTriplet(species1, p3Angle, angle);
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
                pmBonding.setBondingPotentialQuad(species1, p4BondTorsionArray2[i], torsion);
                break;
            }
        }
    }
    public void setBondStretchTraPPE(ISpecies species, List<int[]> bond, PotentialMasterBonding pmBonding, String name){
        double Vi =0, Vj =0, V=0, Vtrue=0,  type;
        int p;
        int i =0;
        int[] bondArr = null;
      //  Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT),Mole.UNIT);
        P2HarmonicUFF p2Bond = new P2HarmonicUFF(Double.POSITIVE_INFINITY,  1.0);
        if(name.equals("ethane") || name.equals("ethene") || name.equals("n2")){
            bondArr = new int[]{0,1};
            bond.add(bondArr);
        } else if (name.equals("propene")){
            bondArr = new int[]{0,1};
            bond.add(bondArr);
            bondArr = new int[]{2,1};
            bond.add(bondArr);
        } else if (name.equals("propane")){
            bondArr = new int[]{0,1};
            bond.add(bondArr);
            bondArr = new int[]{2,1};
            bond.add(bondArr);
        }else if( name.equals("co2")){
            bondArr = new int[]{0,1};
            bond.add(bondArr);
            bondArr = new int[]{2,0};
            bond.add(bondArr);
        } else if (name.equals("ch4")) {
            
        } else {
            System.out.println(name);
            throw new RuntimeException("TraPPE gas unknown");
        }



        //bondParams.add(bondConstant);
        pmBonding.setBondingPotentialPair(species, p2Bond, bond );
    }
    void addUniqueElements(Set<AtomType> set, List<AtomType> list) {
        for (AtomType atomType : list) {
            if (!set.contains(atomType)) {
                set.add(atomType);
            }
        }
    }

    public void doLJ(List<List<AtomType>> pairsAtoms, PotentialMasterCell potentialMasterCell, LJUFF[] p2LJ, int pairAtomSize, double truncatedRadius){
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        int i = 0;
        UFF uff = new UFF();
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        // System.out.println(pairsAtoms + " Pairs");
        P2SoftSphericalSumTruncated[] p2Trunc = new P2SoftSphericalSumTruncated[pairAtomSize];
        for(List<AtomType>individualPair: pairsAtoms){
            AtomType atomNameOne = individualPair.get(0);
            AtomType atomNameTwo = individualPair.get(1);
            String atomTypeStringOne = String.valueOf(atomNameOne);
            String atomTypeStringTwo = String.valueOf(atomNameTwo);
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);
            double[] iKey = pdbReaderMOP.atomicPot(atomTypeOne);
            double[] jKey = pdbReaderMOP.atomicPot(atomTypeTwo);
         /*   if (atomTypeOne.equals("CX") && atomTypeTwo.equals("CX")){
                System.out.println("Here in CXCX");
                continue;
            }*/
            // System.out.println(atomTypeOne+ " "+ Arrays.toString(iKey)+" " +atomTypeTwo+ " "+ Arrays.toString(jKey));
            if(iKey == null){
                iKey = returnGaFF(atomTypeOne);
                // System.out.println("Null IKEY "+Arrays.toString(iKey));
                epsilonIKey = kcals.toSim(iKey[0]);
                sigmaIKey = iKey[1];
            } else {
                //System.out.println(Arrays.toString(iKey));
                epsilonIKey = kcals.toSim(iKey[3]);
                sigmaIKey = iKey[2];
            }
            if(jKey == null){
                jKey = returnGaFF(atomTypeTwo);
                // System.out.println("Null JKEY "+ Arrays.toString(jKey));
                epsilonJKey = kcals.toSim(iKey[0]);
                sigmaJKey = jKey[1];
            } else {
                // System.out.println(Arrays.toString(jKey));
                epsilonJKey = kcals.toSim(jKey[3]);
                sigmaJKey = jKey[2];
            }
           // System.out.println( atomTypeOne+ " "+ atomTypeTwo+ " " + sigmaIKey + " "+ sigmaJKey + " " + epsilonIKey + " " +epsilonJKey);
            p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);
            //System.out.println(" \n");
            p2Trunc[i] = new P2SoftSphericalSumTruncated(truncatedRadius, p2LJ[i]);
            potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2Trunc[i], new double[]{1, 0, 0, 1});
            //potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2LJ[i], new double[]{1, 0, 0, 1}, truncatedRadius);
            i++;
        }
    }
   /* public void doLJ(List<List<AtomType>> pairsAtoms, PotentialMasterCell potentialMasterCell, LJUFF[] p2LJ, int pairAtomSize, double truncatedRadius, boolean ifTraPPE){
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        SpeciesGasTraPPE speciesTraPPE = new SpeciesGasTraPPE();
        int i = 0;
        UFF uff = new UFF();
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        // System.out.println(pairsAtoms + " Pairs");
        P2SoftSphericalSumTruncated[] p2Trunc = new P2SoftSphericalSumTruncated[pairAtomSize];
        for(List<AtomType>individualPair: pairsAtoms){
            AtomType atomNameOne = individualPair.get(0);
            AtomType atomNameTwo = individualPair.get(1);
            String atomTypeStringOne = String.valueOf(atomNameOne);
            String atomTypeStringTwo = String.valueOf(atomNameTwo);
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);
            double[] iKey = pdbReaderMOP.atomicPot(atomTypeOne);
            double[] jKey = pdbReaderMOP.atomicPot(atomTypeTwo);
            if (atomTypeOne.equals("CX") && atomTypeTwo.equals("CX")){
                System.out.println("Here in CXCX");
                continue;
            }
            // System.out.println(atomTypeOne+ " "+ Arrays.toString(iKey)+" " +atomTypeTwo+ " "+ Arrays.toString(jKey));
            if(iKey == null){
                iKey = returnGaFF(atomTypeOne);
                // System.out.println("Null IKEY "+Arrays.toString(iKey));
                epsilonIKey = kcals.toSim(iKey[1]);
                sigmaIKey = iKey[0];
            } else if (iKey == null && ifTraPPE) {
                iKey  = speciesTraPPE.atomicPot(atomTypeOne);
                sigmaIKey = iKey[0];
                epsilonIKey = Math.pow(2,1.0/6.0)*0.001987*Kelvin.UNIT.toSim(iKey[1]);
            } else {
                //System.out.println(Arrays.toString(iKey));
                epsilonIKey = kcals.toSim(iKey[3]);
                sigmaIKey = iKey[2];
            }
            if(jKey == null && !ifTraPPE){
                jKey = returnGaFF(atomTypeTwo);
                // System.out.println("Null JKEY "+ Arrays.toString(jKey));
                epsilonJKey = kcals.toSim(iKey[1]);
                sigmaJKey = jKey[0];
            }else if (ifTraPPE) {
                jKey  = speciesTraPPE.atomicPot(atomTypeTwo);
                sigmaJKey = jKey[0];
                epsilonJKey =Math.pow(2,1.0/6.0)*0.001987*Kelvin.UNIT.toSim(jKey[1]);
            } else {
                // System.out.println(Arrays.toString(jKey));
                epsilonJKey = kcals.toSim(jKey[3]);
                sigmaJKey = jKey[2];
            }
            // System.out.println( atomTypeOne+ " "+ atomTypeTwo+ " " + sigmaIKey + " "+ sigmaJKey + " " + epsilonIKey + " " +epsilonJKey);
            p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);
            //System.out.println(" \n");
            p2Trunc[i] = new P2SoftSphericalSumTruncated(truncatedRadius, p2LJ[i]);
            potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2Trunc[i], new double[]{1, 0, 0, 1});
            //potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2LJ[i], new double[]{1, 0, 0, 1}, truncatedRadius);
            i++;
        }
    }*/

    public void doLJElectrostatic(List<List<AtomType>> pairsAtoms, PotentialMasterCell potentialMasterCell, LJUFF[] p2LJ, P2Electrostatic[] P2Electrostatics, IPotential2[] p2lj, int pairAtomSize, double truncatedRadius, boolean doElectrostatics){
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        int i = 0;
        UFF uff = new UFF();
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        // System.out.println(pairsAtoms + " Pairs");
        P2SoftSphericalSumTruncated[] p2Trunc = new P2SoftSphericalSumTruncated[pairAtomSize];
        for(List<AtomType>individualPair: pairsAtoms){
            AtomType atomNameOne = individualPair.get(0);
            AtomType atomNameTwo = individualPair.get(1);
            String atomTypeStringOne = String.valueOf(atomNameOne);
            String atomTypeStringTwo = String.valueOf(atomNameTwo);
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);
           // System.out.println(atomTypeOne+ " " +atomTypeTwo);
            double[] iKey = pdbReaderMOP.atomicPot(atomTypeOne);
            double[] jKey = pdbReaderMOP.atomicPot(atomTypeTwo);

            if(iKey == null){
                iKey = returnGaFF(atomTypeOne);
                // System.out.println("Null IKEY "+Arrays.toString(iKey));
                epsilonIKey = kcals.toSim(iKey[1]);
                sigmaIKey = iKey[0];
            } else {
                //System.out.println(Arrays.toString(iKey));
                epsilonIKey = kcals.toSim(iKey[3]);
                sigmaIKey = iKey[2];
            }
            if(jKey == null){
                jKey = returnGaFF(atomTypeTwo);
                // System.out.println("Null JKEY "+ Arrays.toString(jKey));
                epsilonJKey = kcals.toSim(iKey[0]);
                sigmaJKey = jKey[1];
            } else {
                // System.out.println(Arrays.toString(jKey));
                epsilonJKey = kcals.toSim(jKey[3]);
                sigmaJKey = jKey[2];
            }
            if(doElectrostatics){
                double chargeOne = pdbReaderMOP.getatomCharge(atomTypeOne);
                double chargeTwo = pdbReaderMOP.getatomCharge(atomTypeTwo);
                P2Electrostatics[i] = uff.electroUFF(chargeOne, chargeTwo);
            }

            // System.out.println( atomTypeOne+ " "+ atomTypeTwo+ " " + sigmaIKey + " "+ sigmaJKey + " " + epsilonIKey + " " +epsilonJKey);
            p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);

            if(doElectrostatics){
                p2Trunc[i] = new P2SoftSphericalSumTruncatedForceShifted(truncatedRadius, p2LJ[i], P2Electrostatics[i]);
            }else {
                TruncationFactory tf = new TruncationFactoryForceShift(truncatedRadius);
                p2lj[i] = tf.make(p2LJ[i]);
            }
            //System.out.println(" \n");
            potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2lj[i], new double[]{1, 0, 0, 1});
            //potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2LJ[i], new double[]{1, 0, 0, 1}, truncatedRadius);
            i++;
        }
    }
    public void doLJElectrostatic(List<List<AtomType>> pairsAtoms, PotentialMasterCell potentialMasterCell, LJUFF[] p2LJ, P2Electrostatic[] P2Electrostatics, int pairAtomSize, double truncatedRadius, boolean doElectrostatics, boolean ifTraPPE){
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        SpeciesGasTraPPE speciesTraPPE = new SpeciesGasTraPPE();
        int i = 0;
        UFF uff = new UFF();
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        // System.out.println(pairsAtoms + " Pairs");
        P2SoftSphericalSumTruncated[] p2Trunc = new P2SoftSphericalSumTruncated[pairAtomSize];
        for(List<AtomType>individualPair: pairsAtoms){
            AtomType atomNameOne = individualPair.get(0);
            AtomType atomNameTwo = individualPair.get(1);
            String atomTypeStringOne = String.valueOf(atomNameOne);
            String atomTypeStringTwo = String.valueOf(atomNameTwo);
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);

            double[] iKey = pdbReaderMOP.atomicPot(atomTypeOne);
            double[] jKey = pdbReaderMOP.atomicPot(atomTypeTwo);
           /* if (atomTypeOne.equals("CX") && atomTypeTwo.equals("CX")){
                System.out.println("Here in CXCX");
                continue;
            }*/
            // System.out.println(atomTypeOne+ " "+ Arrays.toString(iKey)+" " +atomTypeTwo+ " "+ Arrays.toString(jKey));
            if(iKey == null && ifTraPPE){
                iKey  = speciesTraPPE.atomicPot(atomTypeOne);
                sigmaIKey = iKey[0];
              //  double val = Math.pow(2,1.0/6.0);
              //  epsilonIKey = Math.pow(2,1.0/6.0)*0.001987*Kelvin.UNIT.toSim(iKey[1]);
               // epsilonIKey = Math.pow(2,1.0/6.0)*Kelvin.UNIT.toSim(iKey[1]);
                epsilonIKey = Kelvin.UNIT.toSim(Math.pow(2,1.0/6.0)*iKey[1]);
            } else if (iKey == null ) {
                iKey = returnGaFF(atomTypeOne);
                // System.out.println("Null IKEY "+Arrays.toString(iKey));
                epsilonIKey = kcals.toSim(iKey[0]);
                sigmaIKey = iKey[1];
            } else {
                //System.out.println(Arrays.toString(iKey));
                epsilonIKey = kcals.toSim(iKey[3]);
                sigmaIKey = iKey[2];
            }
            if(jKey == null && !ifTraPPE){
                jKey = returnGaFF(atomTypeTwo);
                // System.out.println("Null JKEY "+ Arrays.toString(jKey));
                epsilonJKey = kcals.toSim(iKey[0]);
                sigmaJKey = jKey[1];
            }else if (ifTraPPE) {
                jKey  = speciesTraPPE.atomicPot(atomTypeTwo);
                sigmaJKey = jKey[0];
              //  epsilonJKey =Math.pow(2,1.0/6.0)*0.001987*Kelvin.UNIT.toSim(jKey[1]);
               // epsilonJKey =Math.pow(2,1.0/6.0)*Kelvin.UNIT.toSim(jKey[1]);
                epsilonJKey =Kelvin.UNIT.toSim(Math.pow(2,1.0/6.0)*jKey[1]);
            } else {
                // System.out.println(Arrays.toString(jKey));
                epsilonJKey = kcals.toSim(jKey[3]);
                sigmaJKey = jKey[2];
            }
            if(doElectrostatics){
                double chargeOne, chargeTwo = 0.0;
                if(ifTraPPE){
                    chargeOne = speciesTraPPE.getCharge()[0];
                    chargeTwo = speciesTraPPE.getCharge()[0];
                } else {
                    chargeOne = pdbReaderMOP.getatomCharge(atomTypeOne);
                    chargeTwo = pdbReaderMOP.getatomCharge(atomTypeTwo);
                }

                P2Electrostatics[i] = uff.electroUFF(chargeOne, chargeTwo);
            }

            // System.out.println( atomTypeOne+ " "+ atomTypeTwo+ " " + sigmaIKey + " "+ sigmaJKey + " " + epsilonIKey + " " +epsilonJKey);
            p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);

            if(doElectrostatics){
                p2Trunc[i] = new P2SoftSphericalSumTruncated(truncatedRadius, p2LJ[i], P2Electrostatics[i]);
            }else {
                p2Trunc[i] = new P2SoftSphericalSumTruncated(truncatedRadius, p2LJ[i]);
            }
            //System.out.println(" \n");
            potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2Trunc[i], new double[]{1, 0, 0, 1});
            //potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2LJ[i], new double[]{1, 0, 0, 1}, truncatedRadius);
            i++;
        }
    }

    public void doLJ(List<List<AtomType>> pairsAtoms, LJUFF[] p2LJ, IPotential2[] p2lj, double rc, PotentialMasterCell potentialMasterCell, double[] sigmaIJ){
        PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
        int i = 0;
        UFF uff = new UFF();
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
       // System.out.println(pairsAtoms + " Pairs");
        for(List<AtomType>individualPair: pairsAtoms){
            AtomType atomNameOne = individualPair.get(0);
            AtomType atomNameTwo = individualPair.get(1);
            String atomTypeStringOne = String.valueOf(atomNameOne);
            String atomTypeStringTwo = String.valueOf(atomNameTwo);
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);
            double[] iKey = pdbReaderReplica.atomicPot(atomTypeOne);
            double[] jKey = pdbReaderReplica.atomicPot(atomTypeTwo);
          //  System.out.println(atomTypeOne+" "+Arrays.toString(iKey) + " " +atomTypeTwo+" "+ Arrays.toString(jKey));
            double epsilonIKey = kcals.toSim(iKey[3]);
            double epsilonJKey = kcals.toSim(jKey[3]);
           // System.out.println(iKey[3] + " " +jKey[3]);
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
    public void doLJGrapheneMixed(List<List<AtomType>> pairsAtoms, LJUFF[] p2LJ, IPotential2[] p2lj, P2Electrostatic[] P2Electrostatics, int rc,  PotentialMaster potentialMaster,PotentialMasterCell potentialMasterCell, double[] sigmaIJ, boolean doElectrostatics){
        PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        grapheneReader grapheneReader = new grapheneReader();
        P2SoftSphericalSumTruncated[] p2Trunc = new P2SoftSphericalSumTruncated[pairsAtoms.size()];
        int i = 0;
        UFF uff = new UFF();
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        for(List<AtomType>individualPair: pairsAtoms){
            AtomType atomNameOne = individualPair.get(0);
            AtomType atomNameTwo = individualPair.get(1);
            String atomTypeStringOne = String.valueOf(atomNameOne);
            String atomTypeStringTwo = String.valueOf(atomNameTwo);
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);
            double[] iKey = pdbReaderReplica.atomicPot(atomTypeOne);
            double[] jKey = pdbReaderReplica.atomicPot(atomTypeTwo);
            System.out.println(atomTypeOne+" "+Arrays.toString(iKey) + " " +atomTypeTwo+" "+ Arrays.toString(jKey));
            double epsilonIKey = kcals.toSim(iKey[3]);
            double epsilonJKey = kcals.toSim(jKey[3]);
            double sigmaIKey = iKey[2];
            double sigmaJKey = jKey[2];
            if(doElectrostatics){
               // System.out.println(atomTypeOne + " " + atomTypeTwo);
                double chargeOne = pdbReaderMOP.getatomCharge(atomTypeOne);
                double chargeTwo = pdbReaderMOP.getatomCharge(atomTypeTwo);
                P2Electrostatics[i] = uff.electroUFF(chargeOne, chargeTwo);
              //  System.out.println(atomTypeOne + "  " + chargeOne + " " + atomTypeTwo + " " + chargeTwo);
            }

            sigmaIJ[i] = (sigmaIKey + sigmaJKey) / 2;
            TruncationFactory tf = new TruncationFactoryForceShift(rc);
            p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);
            if(doElectrostatics){
                p2Trunc[i] = new P2SoftSphericalSumTruncated(rc, p2LJ[i], P2Electrostatics[i]);
            }else {
                p2Trunc[i] = new P2SoftSphericalSumTruncated(rc, p2LJ[i]);
            }
            potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2Trunc[i], new double[]{1, 0, 0, 1});
            potentialMaster.setPairPotential(atomNameOne, atomNameTwo, p2Trunc[i], new double[]{1, 0, 0, 1});
            //potentialMaster.setPairPotential(atomNameOne, atomNameTwo, p2lj[i], new double[]{1, 0, 0, 1});
            //potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2lj[i], new double[]{1, 0, 0, 1});
            i++;
        }
    }

    public void doLJMD(List<List<AtomType>> pairsAtoms, LJUFF[] p2LJ, IPotential2[] p2lj, int rc, PotentialMaster potentialMaster, double[] sigmaIJ, PotentialMasterCell potentialMasterCell){
        PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
        int i = 0;
        UFF uff = new UFF();
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
       // System.out.println(pairsAtoms + " Pairs");
        for(List<AtomType>individualPair: pairsAtoms){
            AtomType atomNameOne = individualPair.get(0);
            AtomType atomNameTwo = individualPair.get(1);
            String atomTypeStringOne = String.valueOf(atomNameOne);
            String atomTypeStringTwo = String.valueOf(atomNameTwo);
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);
            //System.out.println(atomTypeOne + " " + atomTypeTwo);
            double[] iKey = pdbReaderReplica.atomicPot(atomTypeOne);
            double[] jKey = pdbReaderReplica.atomicPot(atomTypeTwo);
            double epsilonIKey = kcals.toSim(iKey[3]);
            double epsilonJKey = kcals.toSim(jKey[3]);
            double sigmaIKey = iKey[2];
            double sigmaJKey = jKey[2];
            sigmaIJ[i] = (sigmaIKey + sigmaJKey) / 2;
            TruncationFactory tf = new TruncationFactoryForceShift(rc);
            p2LJ[i] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);
            p2lj[i] = tf.make(p2LJ[i]);
            // potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo,p2lj[i]);
            potentialMaster.setPairPotential(atomNameOne, atomNameTwo, p2lj[i], new double[]{1, 0, 0, 1});
            potentialMasterCell.setPairPotential(atomNameOne, atomNameTwo, p2lj[i], new double[]{1, 0, 0, 1});
            i++;
        }
    }

    public double[] returnGaFF (String atomTypeOne){
        grapheneReader grapheneReader = new grapheneReader();
        double[] num = grapheneReader.atomicPotGraphene(atomTypeOne);
        if (num == null){
            throw new RuntimeException("atomType One " +atomTypeOne);

        }
        return num;
    }
    public List<List<AtomType>> getSpeciesPairs (ISpecies species){
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
    public List<List<AtomType>> listFinal(List<AtomType> list1){
        List<List<AtomType>> pairsAtomsTotal = new ArrayList<>();
        for(i=0; i<list1.size(); i++) {
            for (int j = 0; j < list1.size(); j++) {
                if(i<=j){
                    List<AtomType> subPair = new ArrayList<>();
                    subPair.add(list1.get(i));
                    subPair.add(list1.get(j));
                    pairsAtomsTotal.add(subPair);
                }
            }
        }
        return pairsAtomsTotal;
    }
    public List<List<AtomType>> listFinal(Set<AtomType> set) {
        List<AtomType> list1 = new ArrayList<>(set);
        List<List<AtomType>> pairsAtomsTotal = new ArrayList<>();
        for (int i = 0; i < list1.size(); i++) {
            for (int j = 0; j < list1.size(); j++) {
                if (i <= j) {
                    List<AtomType> subPair = new ArrayList<>();
                    subPair.add(list1.get(i));
                    subPair.add(list1.get(j));
                    pairsAtomsTotal.add(subPair);
                }
            }
        }
        return pairsAtomsTotal;
    }

    public List<List<AtomType>> listFinal(Set<AtomType> set, Set<AtomType> membrane) {
        List<AtomType> list1 = new ArrayList<>(set);
        List<AtomType> list2 = new ArrayList<>(membrane);
        List<List<AtomType>> pairsAtomsTotal = new ArrayList<>();
        for (int i = 0; i < list2.size(); i++) {
            for (int j = 0; j < list1.size(); j++) {
                List<AtomType> subPair = new ArrayList<>();
                subPair.add(list2.get(i));
                subPair.add(list1.get(j));
                pairsAtomsTotal.add(subPair);
            }
        }
        return pairsAtomsTotal;
    }
    public List<AtomType> listTwoSpeciesPairs(List<AtomType> list1, List<AtomType> list2){
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
        return list3;
    }

    public List<List<AtomType>> listGrapheneSpecial (List<AtomType> listGraph, List<AtomType> listGas){
        List<List<AtomType>> pairsAtomsTotal = new ArrayList<>();
        for(int i=0; i<listGraph.size(); i++){
            for(int j=0; j<listGas.size(); j++){
                List<AtomType> subPair = new ArrayList<>();
                subPair.add(listGraph.get(i));
                subPair.add(listGas.get(j));
                pairsAtomsTotal.add(subPair);
            }
        }
        //System.out.println(pairsAtomsTotal);
        return pairsAtomsTotal;
    }
    public IPotential2[][] makeAtomPotentials(SpeciesManager sm) {
        // we could try to store the potentials more compactly, but it doesn't really matter
        ISpecies species = sm.getSpecies(sm.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
       // System.out.println(lastTypeIndex + 1+ " "+lastTypeIndex + 1 + " lastTypeIndex" + " "+species.getAtomType(species.getUniqueAtomTypeCount() - 1));
        return new IPotential2[lastTypeIndex + 1][lastTypeIndex + 1];
    }


    public void setupBondStrechOPLS(ISpecies species){}

    public void setupBondStrech(ISpecies species1,Map<String[],List<int[]>> bondTypesMap1, Map<String[],List<int[]>> angleTypesMap1,Map<String[],List<int[]>> torsionTypesMap1,ArrayList<Integer> bondsNum1,ArrayList<Integer> bondList1, List<int[]>quadrupletsSorted1, Map<Integer, String> atomIdentifierMapModified1,Map<String, double[]> atomicPotMap1, PotentialMasterBonding.FullBondingInfo bondingInfo1, PotentialMasterBonding pmBonding  ){
        double Vi =0, Vj =0, V=0, Vtrue=0,  type;
        int p;
        int i =0;
        double bondOrder;
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        for (Map.Entry<String[], List<int[]>> entry : bondTypesMap1.entrySet()) {
            String[] bondType = entry.getKey();
            List<int[]> bonds = entry.getValue();
            // System.out.println(Arrays.toString(bondType) + ": " + Arrays.deepToString(bonds.toArray()));
            for(int[]bondIndividual: bonds){
             //  bonds.add(bondIndividual);
                double[] bondParamsArray = new double[2];
                double[] bondConstant = new double[2];
                atom1 = bondIndividual[0];
                atom2 = bondIndividual[1];
                atomName1 = atomIdentifierMapModified1.get(atom1);
                atomName2 = atomIdentifierMapModified1.get(atom2);
                double[] atomOnePot = atomicPotMap1.get(atomName1);
                double[] atomTwoPot = atomicPotMap1.get(atomName2);
                if(atomName1.equals("C_2p") && atomName2.equals("H_p")){
                     bondOrder =1;
                }else if (atomName1.equals("C_2p") && atomName2.equals("C_2p")){
                    bondOrder =2;
                } else if (atomName1.equals("C_3p") && atomName2.equals("C_2p")) {
                    bondOrder=1;
                }else if (atomName1.equals("C_3") && atomName2.equals("O_1")) {
                    bondOrder=2;
                } else {
                    bondOrder = bondsNum1.get(i);
                }
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
               // System.out.println(atomName1 + " " +atomName2 + " " +bondOrder);

              /*  if(atomName1.equals("C_Ar") && atomName2.equals("C_Ar")){
                    bondOrder = 1.5;
                } else {
                    bondOrder = 1;
                }*/
              //  System.out.println(bondOrder +" bondorder");
                //  System.out.println(Arrays.toString(dupletsSorted.get(i)) + " " + bondOrder+ " " + atomName1 + " " + atomName2+" "+ Arrays.toString(atomOnePot) +" " + Arrays.toString(atomTwoPot));
                bondParamsArray= UFF.bondUFF (atomOnePot[0],  atomTwoPot[0], atomOnePot[5],  atomTwoPot[5], atomOnePot[6], atomTwoPot[6], bondOrder);
                //System.out.println(Arrays.toString(bondParamsArray) + " ArrayToString");
                bondConstant = UFF.BondConstantArray(bondParamsArray[0], bondParamsArray[1]);
                //System.out.println(Arrays.toString(bondConstant) + " arrayConstant");
                P2HarmonicUFF p2Bond = new P2HarmonicUFF(bondParamsArray[0],  bondParamsArray[1]);
                //bondParams.add(bondConstant);
                pmBonding.setBondingPotentialPair(species1, p2Bond, bonds);
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
                P3BondAngleUFF p3Angle = new P3BondAngleUFF(angleParamsArray[0],  angleParamsArray[1], angleParamsArray[2], angleParamsArray[3],Degree.UNIT.toSim( atomTwoPot[1]), 0, caseNum);
                pmBonding.setBondingPotentialTriplet(species1, p3Angle, angle);
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
                bondOrder = 1 ;
                torsionParamsArray = UFF.torsionUFF(Vtrue, p, bondOrder);
                p4BondTorsionArray2[i] = new P4BondTorsionUFF(torsionParamsArray[0], (int) torsionParamsArray[1], torsionParamsArray[2]);
                double[] array = {torsionParamsArray[0], torsionParamsArray[1], torsionParamsArray[2]};
                p4ValueArray2.add(array);
                //System.out.println(Arrays.toString(array));
                pmBonding.setBondingPotentialQuad(species1, p4BondTorsionArray2[i], torsion);
                i++;
            }
        }
    }

    public void setupBondStretchOPLS(ISpecies species1,Map<String[],List<int[]>> bondTypesMap1, Map<String[],List<int[]>> angleTypesMap1,Map<String[],List<int[]>> torsionTypesMap1,ArrayList<Integer> bondsNum1,ArrayList<Integer> bondList1, List<int[]>quadrupletsSorted1, Map<Integer, String> atomIdentifierMapModified1,Map<String, double[]> atomicPotMap1, PotentialMasterBonding.FullBondingInfo bondingInfo1, PotentialMasterBonding pmBonding  ){
        double Vi =0, Vj =0, V=0, Vtrue=0,  type;
        int p;
        int i =0;
        double bondOrder;
        PDBDataExtracterOPLS pdbDataExtracterOPLS = new PDBDataExtracterOPLS();
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        for (Map.Entry<String[], List<int[]>> entry : bondTypesMap1.entrySet()){
            String[] bondType = entry.getKey();
            List<int[]> bonds = entry.getValue();
            for(int[]bondIndividual: bonds){
                double[] bondParamsArray = new double[2];
                double[] bondConstant = new double[2];
                atom1 = bondIndividual[0];
                atom2 = bondIndividual[1];
                atomName1 = atomIdentifierMapModified1.get(atom1);
                atomName2 = atomIdentifierMapModified1.get(atom2);
                bondParamsArray = pdbDataExtracterOPLS.bondValueSender(atomName1, atomName2);
            }
        }
    }
}
