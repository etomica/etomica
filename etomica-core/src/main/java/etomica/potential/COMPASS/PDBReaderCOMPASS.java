package etomica.potential.COMPASS;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.chem.elements.*;
import etomica.molecule.IMolecule;
//import etomica.potential.GAFF.PDBReaderGAFF;
import etomica.potential.UFF.PDBReaderMOP;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;

import java.util.*;
/*
public class PDBReaderCOMPASS {
    ISpecies species;
    Map<Integer, Integer> assignCharges = new HashMap<>();
    Map<Integer, AtomType> atomIdentifierCOMPASS = new HashMap<>();
    static Map<Integer, String> modifiedAtomIdentifierMap = new HashMap<>();
    public static Map<String, AtomType> elementReceiverMap = new HashMap<>();
    public static Map<Integer, Vector> positions = new HashMap<>();
    public Map<Double, Double> chargeMap = new HashMap<>();
    public Map<String, AtomType> typeMap = new HashMap<>();
    public Map<String, AtomType> typeMapNew = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectedAtoms = new ArrayList<>();
    public Map<Integer,String> atomIdentifierMapMod= new HashMap<>();
    public ISpecies getSpecies (String confName, boolean isDynamic, boolean setMOPmassInfinite, Vector centreMOP){
        double massSum = 0;
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        SpeciesBuilder speciesBuilderNewMod =  new SpeciesBuilder(Space3D.getInstance());
        AtomType typeNew;
        PDBReaderMOP PDBReader = new PDBReaderMOP();
        PDBReaderCOMPASS PDBReaderCOMPASS = new PDBReaderCOMPASS();
        PDBReader.readPDBFile(confName);
        ArrayList<ArrayList<Integer>> connectivity = PDBReader.getConnectivity();
        Map<Integer, String> atomMap = PDBReader.getAtomMap();
        ArrayList<Integer> bondList =PDBReader.setBondList(connectivity, atomMap);
        //System.out.println(bondList + " bondList");
        System.out.println(atomMap+ "  " + atomMap.size());
        atomIdentifierMapMod = atomIdentifierMapModifiedCOMPASS(connectivity, atomMap);
        // System.out.println(atomIdentifierMapModified + "atomIdentifierMapModified");
         connectedAtoms =PDBReader.getConnectivity();
        //System.out.println(connectedAtoms+ ": connectedAtom");
        ArrayList<ArrayList<Integer>> connectivityModified = PDBReader.getconnectivityModified(connectedAtoms);
        // System.out.println(connectivityModified+ ": connectedAtomModified" );
        // System.out.println(atomMap + ": atomMap");
        Map<Integer, String> atomMapModified =  PDBReader.getatomMapModified(atomMap);
        // System.out.println(atomMapModified + ": atomMapModified");
        positions = PDBReader.getPositions();
        System.out.println(positions);
        System.out.println(connectedAtoms +" O");
        System.out.println(connectivityModified+" O");
        System.out.println(atomMapModified+" O");
        System.out.println(atomIdentifierMapMod+" O");
        ArrayList<Integer> bondsNum = PDBReader.bondsAmongAtoms();
        System.out.println(bondsNum + " bond Amongatoms");
        List<int[]> duplets =  PDBReader.getBondedAtomList(connectivityModified);

        List<int[]> listOfAngleModified =  PDBReader.getAngleList(connectivityModified);
        List<String[]> angleSorted =  PDBReaderCOMPASS.getAngleSorted(listOfAngleModified, modifiedAtomIdentifierMap);
        System.out.println(Arrays.deepToString(angleSorted.toArray()));
       // List<int[]> dupletsSorted = PDBReader.bondSorter(duplets, atomIdentifierMapMod);
        //System.out.println(Arrays.deepToString(dupletsSorted.toArray()) + " duplets Sorted");
        //System.out.println(Arrays.deepToString(duplets.toArray())+ ": listOfBonds");
       // List<int[]> tripletsSorted = PDBReader.angleSorter(listOfAngleModified, atomIdentifierMapMod);
        List<int[]> listOfTorsionModified =  PDBReader.getTorsionList(connectivity);
        //System.out.println( Arrays.deepToString(listOfTorsionModified.toArray()) + "torsionModified");
      //  List<int[]> quadrupletsSorted = PDBReader.torsionSorter(listOfTorsionModified, atomIdentifierMapMod);
        // System.out.println(Arrays.deepToString(quadrupletsSorted.toArray()) + "quadrupletsSorted");
       // Map<String[],List<int[]>> torsionTypesMap=  PDBReader.idenTorsionTypes(quadrupletsSorted, atomIdentifierMapMod);
        System.out.println(Arrays.deepToString(listOfAngleModified.toArray()));
        List<int[]> bondBonds = PDBReaderCOMPASS.getBondBond(listOfAngleModified);
       /* List<int[]> bondAltBond = PDBReaderCOMPASS.getBondAltBond(listOfTorsionModified);
        List<int[]> bondAngle = PDBReaderCOMPASS.getBondAngle(listOfAngleModified);
        List<int[]> angleAngle = PDBReaderCOMPASS.getAngleAngle(listOfTorsionModified);
        List<int[]> endBondTorsion = PDBReaderCOMPASS.getEndBondTorsion(listOfTorsionModified);
        List<int[]> midBondTorsion = PDBReaderCOMPASS.getMidBondTorsion(listOfTorsionModified);
        List<int[]> angleTorsion = PDBReaderCOMPASS.getAngleTorsion(listOfTorsionModified);
        List<int[]> angleAngleTorsion = PDBReaderCOMPASS.getAngleAngleTorsion(listOfTorsionModified);
        System.out.println(Arrays.deepToString(listOfInversions.toArray()) + " in Main" );
        System.out.println( PDBReader.positions + " Here positions");
        Space space = Space3D.getInstance();
        Vector center = space.makeVector();
        Vector dr = Vector.d(center.getD());
        System.out.println(atomIdentifierMapMod);
        setConnectivity(connectedAtoms);
        setAtomIdentifierModified(atomIdentifierMapMod);
        //   System.out.println(atomIdentifierMap);
        for(int i = 0; i < atomIdentifierMapMod.size(); i++) {
            String nameNew = String.valueOf(atomIdentifierMapMod.get(i));
            // System.out.println(atomIdentifierMap.get(i));
            AtomType newName = returnElementCOMPASS(nameNew, setMOPmassInfinite);
            // AtomType newNameNew = returnElement(nameNew);
            if (typeMapNew.containsKey(nameNew)) {
                typeNew = typeMapNew.get(nameNew);
            } else {
                typeNew = newName;
                typeMapNew.put(nameNew, typeNew);
            }
            Vector position = positions.get(i+1);
            speciesBuilderNew.addAtom(typeNew, position,  "");
        }
        species= speciesBuilderNew.setDynamic(isDynamic).build();
        if(!setMOPmassInfinite){
            //   System.out.println(typeMapNew + " typeMapNew");

            //   System.out.println(species.getMass() + " first");
            IMolecule molecule = species.makeMolecule();
            IAtomList children = molecule.getChildList();
            int nAtoms = children.size();
            for (int i = 0; i < nAtoms; i++) {
                IAtom a = children.get(i);
                // System.out.println(a.getPosition() + " "+ i);
                double mass = a.getType().getMass();
                if (massSum == 0) {
                    centreMOP.PEa1Tv1(mass, a.getPosition());
                } else {
                    // sum = sum + mass*((sum/n)+pbc(r - sum/n))
                    dr.E(a.getPosition());
                    centreMOP.PEa1Tv1(mass, dr);
                }
                massSum += mass;
            }
            centreMOP.TE(1.0 / massSum);
            System.out.println(centreMOP + " 1 ");
            System.out.println("Part 1");
            for(int i=0; i<atomIdentifierMapMod.size(); i++){
                IAtom a = children.get(i);
                String nameNew = String.valueOf(atomIdentifierMapMod.get(i));
                typeNew = typeMapNew.get(nameNew);
                Vector v = a.getPosition();
                v.ME(centreMOP);
                // System.out.println(v);
                speciesBuilderNewMod.addAtom(typeNew,v, "" );
            }
            species = speciesBuilderNewMod.setDynamic(isDynamic).build();
            System.out.println(species.getMass() + " Second");
        }
        return species;
    }

    public void setAtomIdentifierModified(Map<Integer, String> atomIdentifierMapMod) {
        this.atomIdentifierMapMod = atomIdentifierMapMod;
    }

    public void setConnectivity(ArrayList<ArrayList<Integer>> connectivity){
        this.connectedAtoms = connectivity;
    }

    public ArrayList<ArrayList<Integer>> getConnectivity(){
        return connectedAtoms;
    }
    public Map<Integer, String> getAtomIdentifierMapMod(){return atomIdentifierMapMod;}

    private List<String[]> getAngleSorted(List<int[]> listOfAngleModified, Map<Integer, String> modifiedAtomIdentifierMap) {
        Set<String> uniqueArrays = new HashSet<>();
        List<String[]> result = new ArrayList<>();

        for (int[] mapping : listOfAngleModified) {
            if (mapping[0] < mapping[2]) {
                String[] mappedArray = new String[mapping.length];
                for (int i = 0; i < mapping.length; i++) {
                    mappedArray[i] = modifiedAtomIdentifierMap.get(mapping[i]);
                }
                String arrayString = Arrays.toString(mappedArray);
                if (!uniqueArrays.contains(arrayString)) {
                    uniqueArrays.add(arrayString);
                    result.add(mappedArray);
                }
            }
        }

        return result;
    }

    private List<int[]> getBondBond(List<int[]> listOfAngleModified) {
        List<int[]> bondBond = new ArrayList<>();
        int[] num = new int[2];
        int[] numTwo = new int[2];
        String atomName;
        List<String[]> result = new ArrayList<>();

        for (int[] mapping : listOfAngleModified) {
            String[] mappedArray = new String[mapping.length];
            for (int i = 0; i < mapping.length; i++) {
                mappedArray[i] = modifiedAtomIdentifierMap.get(mapping[i]);
            }
            result.add(mappedArray);
        }
        Set<String> uniqueArrays = new HashSet<>();
        List<String[]> resultNew = new ArrayList<>();
        for (int[] mapping : listOfAngleModified) {
            String[] mappedArray = new String[mapping.length];
            for (int i = 0; i < mapping.length; i++) {
                mappedArray[i] = modifiedAtomIdentifierMap.get(mapping[i]);
            }
            String arrayString = Arrays.toString(mappedArray);
            if (!uniqueArrays.contains(arrayString)) {
                uniqueArrays.add(arrayString);
                resultNew.add(mappedArray);
            }
        }
        System.out.println(Arrays.deepToString(resultNew.toArray()));
        System.out.println(Arrays.deepToString(result.toArray()));
        for (int i=0; i<listOfAngleModified.size(); i++){
            boolean isPresent = false;
            boolean isPresentTwo = false;
            int[] angle = listOfAngleModified.get(i);
            num = new int[]{angle[0], angle[1]};
            numTwo = new int[]{angle[1], angle[2]};
            if (num[0] > num[1]) {
                // Swap elements if not
                int temp = num[0];
                num[0] = num[1];
                num[1] = temp;
            }
            if (numTwo[0] > numTwo[1]) {
                // Swap elements if not
                int temp = numTwo[0];
                numTwo[0] = numTwo[1];
                numTwo[1] = temp;
            }
            for (int[] arr : bondBond) {
                if (Arrays.equals(arr, num)) {
                    isPresent = true;
                    break;
                }
            }
            for (int[] arr : bondBond) {
                if (Arrays.equals(arr, numTwo)) {
                    isPresentTwo = true;
                    break;
                }
            }
            if (!isPresent) {
                bondBond.add(num);
                //System.out.println("Added num to listOfNums");
            } else {
                //System.out.println("num already present in listOfNums");
            }
            if (!isPresentTwo) {
                bondBond.add(numTwo);
               // System.out.println("Added num to listOfNums");
            } else {
               // System.out.println("num already present in listOfNums");
            }

        }
        System.out.println(Arrays.deepToString(bondBond.toArray()));
        return bondBond;
    }

    public static Map<Integer,String> atomIdentifierMapModifiedCOMPASS(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        PDBReaderCOMPASS pdbReaderCOMPASS = new PDBReaderCOMPASS();
        Map<Integer, AtomType> atomIdentifierMap = pdbReaderCOMPASS.atomIdentifierCOMPASS(connectivityModified, atomMapModified);
        for (Map.Entry<Integer, AtomType> entry : atomIdentifierMap.entrySet()) {
            String value = entry.getValue().toString();

            value = value.replace("AtomType[", "").replace("]", "");
            modifiedAtomIdentifierMap.put(entry.getKey(), value);
        }
        System.out.println(modifiedAtomIdentifierMap + " modified Identifier");
        return modifiedAtomIdentifierMap;
    }

    protected Map<Integer, AtomType> atomIdentifierCOMPASS(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified) {
        PDBReaderGAFF readerGAFF = new PDBReaderGAFF();
        if(connectivityModified.size() == 0){
            String element = atomMapModified.get(1);
            if(element.equals("AR")){
                AtomType Ar = new AtomType(Argon.INSTANCE, "Ar");
                atomIdentifierCOMPASS.put(0, Ar);
            } else if (element.equals("NE")) {
                AtomType Ne = new AtomType(Neon.INSTANCE, "Ne");
                atomIdentifierCOMPASS.put(0, Ne);
            } else if (element.equals("HE"))  {
                AtomType He = new AtomType(Helium.INSTANCE, "He");
                atomIdentifierCOMPASS.put(0, He);
            }else if (element.equals("KR"))  {
                AtomType KR = new AtomType(Krypton.INSTANCE, "Kr");
                atomIdentifierCOMPASS.put(0, KR);
            }else if (element.equals("XE"))  {
                AtomType XE = new AtomType(Xenon.INSTANCE, "Xe");
                atomIdentifierCOMPASS.put(0, XE);
            }
        }
        for(int i =0; i<connectivityModified.size(); i++){
            Integer retriveArrayFirstElementNumber = connectivityModified.get(i).get(0);
            String retriveArrayFirstElementName = atomMapModified.get(retriveArrayFirstElementNumber);
            int arrayListSize = connectivityModified.get(i).size();
            int numCarbon=0, numOxy=0, numNitro=0, numHydro=0;
            ArrayList<Integer> innerList=connectivityModified.get(i);
            ArrayList<Integer> innerAtoms= readerGAFF.numAtomCounter(connectivityModified, atomMapModified, innerList);
            numOxy=innerAtoms.get(2);
            numNitro = innerAtoms.get(3);
            numCarbon = innerAtoms.get(0);
            numHydro = innerAtoms.get(1);
            ArrayList<Integer> secondary = readerGAFF.secondaryAtomCounter(connectivityModified,atomMapModified, innerList);
            if(retriveArrayFirstElementName.equals("N")){
                if (arrayListSize == 3 && numNitro == 1) {
                    AtomType N1N = new AtomType(Nitrogen.INSTANCE, "N1N");
                    atomIdentifierCOMPASS.put(i, N1N);
                }else if(numOxy ==1){
                    AtomType N1O = new AtomType(Nitrogen.INSTANCE, "N1O");
                    atomIdentifierCOMPASS.put(i, N1O);
                }else if(numOxy == 2){
                    AtomType N2O = new AtomType(Nitrogen.INSTANCE, "N2O");
                    atomIdentifierCOMPASS.put(i, N2O);
                }
            } else if (retriveArrayFirstElementName.equals("O")) {
                if(arrayListSize == 3 && numNitro == 1 && secondary.get(3) == 1){//NO
                    AtomType O1N = new AtomType(Oxygen.INSTANCE, "O1N");
                    atomIdentifierCOMPASS.put(i, O1N);
                }else if(arrayListSize == 3 && numNitro == 1){//NO2
                    AtomType O1 = new AtomType(Oxygen.INSTANCE, "O1");
                    atomIdentifierCOMPASS.put(i, O1);
                }else if (arrayListSize == 3 && numCarbon == 1 && secondary.get(0) == 1) {//CO
                    AtomType O1C = new AtomType(Oxygen.INSTANCE, "O1C");
                    atomIdentifierCOMPASS.put(i, O1C);
                } else if (arrayListSize == 3 && numCarbon == 1) { //CO2
                    AtomType O2C = new AtomType(Oxygen.INSTANCE, "O2C");
                    atomIdentifierCOMPASS.put(i, O2C);
                } else if (arrayListSize ==3 && numOxy ==1) { //O2
                    AtomType O2O = new AtomType(Oxygen.INSTANCE, "O2O");
                    atomIdentifierCOMPASS.put(i, O2O);
                }else if (arrayListSize ==3 ) { //O=
                    AtomType O2 = new AtomType(Oxygen.INSTANCE, "O2");
                    atomIdentifierCOMPASS.put(i, O2);
                } else if (numCarbon ==2) {//ether
                    AtomType O2E = new AtomType(Oxygen.INSTANCE, "O2E");
                    atomIdentifierCOMPASS.put(i, O2E);
                } else if (numCarbon == 1 && numHydro ==1) {//R-OH
                    AtomType O2H = new AtomType(Oxygen.INSTANCE, "O2H");
                    atomIdentifierCOMPASS.put(i, O2H);
                }
            } else if (retriveArrayFirstElementName.equals("H")) {
                if(numCarbon ==1){
                    AtomType H1 = new AtomType(Hydrogen.INSTANCE, "H1");
                    atomIdentifierCOMPASS.put(i, H1);
                } else if (numHydro ==1) {
                    AtomType H1H = new AtomType(Hydrogen.INSTANCE, "H1H");
                    atomIdentifierCOMPASS.put(i, H1H);
                } else if (numOxy == 1) {
                    AtomType H1O = new AtomType(Hydrogen.INSTANCE, "H1O");
                    atomIdentifierCOMPASS.put(i, H1O);
                }
            } else if (retriveArrayFirstElementName.equals("C")){
                if(arrayListSize == 3 && numOxy ==1){
                    AtomType C1O = new AtomType(Carbon.INSTANCE, "C1O");
                    atomIdentifierCOMPASS.put(i, C1O);
                } else if (arrayListSize == 3 && numOxy ==2) {
                    AtomType C2O = new AtomType(Carbon.INSTANCE, "C2O");
                    atomIdentifierCOMPASS.put(i, C2O);
                }else if (arrayListSize == 4) {
                    AtomType C3A = new AtomType(Carbon.INSTANCE, "C3A");
                    atomIdentifierCOMPASS.put(i, C3A);
                } else if (arrayListSize ==5 && numOxy ==1) {
                    AtomType C4O = new AtomType(Carbon.INSTANCE, "C4O");
                    atomIdentifierCOMPASS.put(i, C4O);
                }else if (arrayListSize == 5){
                    AtomType C4 = new AtomType(Carbon.INSTANCE, "C4");
                    atomIdentifierCOMPASS.put(i, C4);
                }
            }
        }

        return atomIdentifierCOMPASS;
    }

    public AtomType returnElementCOMPASS(String elementName, boolean isInfinite){
        if(isInfinite){
            switch (elementName){
                case "Ar":
                    AtomType Ar = new AtomType(new ElementSimple("Ar", Double.POSITIVE_INFINITY), "Ar");
                    elementReceiverMap.put("Ar", Ar);
                    break;
                case "Cu":
                    AtomType Cu = new AtomType(new ElementSimple("Cu", Double.POSITIVE_INFINITY), "Cu");
                    elementReceiverMap.put("Cu", Cu);
                    break;

            }
        } else {
            elementReceiverMap.put("Ar", new AtomType(Argon.INSTANCE, "Ar"));
            elementReceiverMap.put("HE", new AtomType(Helium.INSTANCE, "HE"));
            elementReceiverMap.put("KR", new AtomType(Krypton.INSTANCE, "KR"));
            elementReceiverMap.put("XE", new AtomType(Xenon.INSTANCE, "XE"));
            elementReceiverMap.put("N1N", new AtomType(Nitrogen.INSTANCE, "N1N"));
            elementReceiverMap.put("N1O", new AtomType(Nitrogen.INSTANCE, "N1O"));
            elementReceiverMap.put("N2O", new AtomType(Nitrogen.INSTANCE, "N2O"));
            elementReceiverMap.put("O1N", new AtomType(Oxygen.INSTANCE, "O1N"));
            elementReceiverMap.put("O1", new AtomType(Oxygen.INSTANCE, "O1"));
            elementReceiverMap.put("O1C", new AtomType(Oxygen.INSTANCE, "O1C"));
            elementReceiverMap.put("O2C", new AtomType(Oxygen.INSTANCE, "O2C"));
            elementReceiverMap.put("O2O", new AtomType(Oxygen.INSTANCE, "O2O"));
            elementReceiverMap.put("O2", new AtomType(Oxygen.INSTANCE, "O2"));
            elementReceiverMap.put("O2E", new AtomType(Oxygen.INSTANCE, "O2E"));
            elementReceiverMap.put("O2H", new AtomType(Oxygen.INSTANCE, "O2H"));
            elementReceiverMap.put("H1", new AtomType(Hydrogen.INSTANCE, "H1"));
            elementReceiverMap.put("H1H", new AtomType(Hydrogen.INSTANCE, "H1H"));
            elementReceiverMap.put("H1O", new AtomType(Hydrogen.INSTANCE, "H1O"));
            elementReceiverMap.put("C1O", new AtomType(Carbon.INSTANCE, "C1O"));
            elementReceiverMap.put("C2O", new AtomType(Carbon.INSTANCE, "C2O"));
            elementReceiverMap.put("C3A", new AtomType(Carbon.INSTANCE, "C3A"));
            elementReceiverMap.put("C4O", new AtomType(Carbon.INSTANCE, "C4O"));
            elementReceiverMap.put("C4", new AtomType(Carbon.INSTANCE, "C4"));
        }
        return elementReceiverMap.get(elementName);
    }

    protected double[] getQuadBondPot(String atomOne, String atomTwo){
        //ArrayList<Integer> quadBondPot = new ArrayList<>();
        double[] quadBondPot = new double[4];
        if(atomOne.equals("C3A") && atomTwo.equals("C3A")){
            quadBondPot = new double[]{1.4170, 470.8361, -627.6179, 1327.6345};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C4")) {
            quadBondPot = new double[]{1.5010,321.9021,-521.8208,572.1628};
        } else if (atomOne.equals("C3A") && atomTwo.equals("H1")) {
            quadBondPot = new double[]{1.0982,372.8251,-803.4526,894.3173};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4")) {
            quadBondPot = new double[]{1.5300,299.6700,-501.7700,679.8100};
        }else if (atomOne.equals("C4") && atomTwo.equals("H1")) {
            quadBondPot = new double[]{1.1010,345.000,-691.8600,844.6000};
        }
        return quadBondPot;
    }

    protected double[] getQuadAnglePot(String atomOne, String atomTwo, String atomThree){
        double[] quadAnglePot = new double[4];
        if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A")){
            quadAnglePot = new double[]{ 118.9000, 61.0226,-34.9931, 0.0000};
        } else if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4")){
            quadAnglePot = new double[]{ 120.0500, 44.7148,-22.7352, 0.0000};
        }else if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("H1")){
            quadAnglePot = new double[]{117.9400, 35.1558,-12.4682, 0.0000};
        }else if(atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("C3A")){
            quadAnglePot = new double[]{ 111.0000, 44.3234,-9.4454, 0.0000};
        }else if(atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("C4")){
            quadAnglePot = new double[]{108.4000, 43.9594,-8.3924,-9.3379};
        }else if(atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("H1")){
            quadAnglePot = new double[]{ 111.0000, 44.3234,-9.4454, 0.0000};
        }else if(atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4")){
            quadAnglePot = new double[]{ 112.6700, 39.5160,-7.4430,-9.5583};
        }else if(atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("H1")){
            quadAnglePot = new double[]{ 110.7700, 41.4530,-10.6040, 5.1290};
        }else if(atomOne.equals("H1") && atomTwo.equals("C4") && atomThree.equals("H1")){
            quadAnglePot = new double[]{ 107.6600, 39.6410,-12.9210,-2.4318};
        }
        return quadAnglePot;
    }

    protected double[] getQuadTorPot(String atomOne, String atomTwo, String atomThree, String atomFour){
        double[] quadTorPot = new double[4];
        if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C3A")){
            quadTorPot = new double[]{ 8.3667, 1.2000, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C4")) {
            quadTorPot = new double[]{ 0.0000, 4.4072, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadTorPot = new double[]{ 0.0000, 3.9661, 0.0000};
        }else if (atomOne.equals("C4") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadTorPot = new double[]{ 0.0000, 1.5590, 0.0000};
        }else if (atomOne.equals("H1") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadTorPot = new double[]{ 0.0000, 2.3500, 0.0000};
        }else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4") && atomFour.equals("C3A")) {
            quadTorPot = new double[]{-0.2802,-0.0678,-0.0122};
        }else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4") && atomFour.equals("C4")) {
            quadTorPot = new double[]{-0.2802,-0.0678,-0.0122};
        }else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadTorPot = new double[]{-0.2802,-0.0678,-0.0122};
        }else if (atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadTorPot = new double[]{-0.0228, 0.0280,-0.1863};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("C4")) {
            quadTorPot = new double[]{ 0.0000, 0.0514,-0.1430};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadTorPot = new double[]{ 0.0000, 0.0316,-0.1681};
        }else if (atomOne.equals("H1") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadTorPot = new double[]{-0.1432, 0.0617,-0.1530};
        }else if ( atomTwo.equals("C3A") && atomThree.equals("C3A")) {
            quadTorPot = new double[]{ 0.0000, 4.5000, 0.0000};
        }else if  (atomTwo.equals("C4") && atomThree.equals("C4")) {
            quadTorPot = new double[]{0.0000, 0.0000,-0.1530};
        }
        return quadTorPot;
    }

    protected double[] getQuadInverPot(String atomOne, String atomTwo, String atomThree, String atomFour){
        //ArrayList<Integer> quadBondPot = new ArrayList<>();
        double[] quadInversPot = new double[4];
        if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C3A")){
            quadInversPot = new double[]{7.1794};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C4")) {
            quadInversPot = new double[]{ 7.8153};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadInversPot = new double[]{4.8912};
        }
        return quadInversPot;
    }

    protected double[] getQuadLJPot(String atomOne){
        //ArrayList<Integer> quadBondPot = new ArrayList<>();
        double[] quadLJPot = new double[2];
        if(atomOne.equals("C3A") ){
            quadLJPot = new double[]{ 3.9150, 0.0680};
        } else if (atomOne.equals("C4") ) {
            quadLJPot = new double[]{3.8540, 0.0620};
        } else if (atomOne.equals("C43") ) {
            quadLJPot = new double[]{3.8540, 0.0400};
        }else if (atomOne.equals("C44") ) {
            quadLJPot = new double[]{ 3.8540, 0.0200};
        }else if (atomOne.equals("H1") ) {
            quadLJPot = new double[]{ 2.8780, 0.0230};
        }
        return quadLJPot;
    }

    protected double[] getBondBondAdjPot(String atomOne, String atomTwo, String atomThree){
        double[] quadBondBondAdj = new double[1];
        if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A")){
            quadBondBondAdj = new double[]{68.2856};
        } else if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4")){
            quadBondBondAdj = new double[]{12.0676};
        }else if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("H1")){
            quadBondBondAdj = new double[]{ 1.0795};
        }else if(atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("H1")){
            quadBondBondAdj = new double[]{ 2.9168};
        }else if(atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("H1")){
            quadBondBondAdj = new double[]{3.3872};
        }else if(atomOne.equals("H1") && atomTwo.equals("C4") && atomThree.equals("H1")){
            quadBondBondAdj = new double[]{5.3316};
        }
        return quadBondBondAdj;
    }

    protected double[] getBondBondAltPot(String atomOne, String atomTwo, String atomThree, String atomFour){
        double[] quadBondBondAdj = new double[1];
        if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C3A")){
            quadBondBondAdj = new double[]{ 53.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{2.5085};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-6.2741};
        }else if (atomOne.equals("C4") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{  0.8743};
        }else if (atomOne.equals("H1") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-1.7077};
        }else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ -3.4826};
        }
        return quadBondBondAdj;
    }

    protected double[] getBondAnglePot(String atomOne, String atomTwo, String atomThree){
        double[] quadBondAngle = new double[2];
        if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A")){
            quadBondAngle = new double[]{ 28.8708, 28.8708};
        } else if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4")){
            quadBondAngle = new double[]{31.0771, 47.0579};
        }else if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("H1")){
            quadBondAngle = new double[]{ 20.0033, 24.2183};
        }else if(atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("H1")){
            quadBondAngle = new double[]{26.4608 ,11.7717};
        }else if(atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4")){
            quadBondAngle = new double[]{ 8.0160 ,8.0160};
        }else if(atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("H1")){
            quadBondAngle = new double[]{ 20.7540, 11.4210};
        }else if(atomOne.equals("H1") && atomTwo.equals("C4") && atomThree.equals("H1")){
            quadBondAngle = new double[]{18.1030, 18.1030};
        }
        return quadBondAngle;
    }

    protected double[] getAngAngPot(String atomOne, String atomTwo, String atomThree, String atomFour){
        double[] quadBondBondAdj = new double[1];
        if(atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C3A") && atomFour.equals("H1")){
            quadBondBondAdj = new double[]{  2.0403};
        } else if (atomOne.equals("H1") && atomTwo.equals("C4") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 3.0118};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-1.8202};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{-0.1729};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-1.3199};
        }else if (atomOne.equals("H1") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-0.4825};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("H1") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{ 1.0827};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("H1") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{2.3794};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("H1") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{0.1184};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("H1") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{0.2738};
        }else if (atomOne.equals("H1") && atomTwo.equals("C4") && atomThree.equals("H1") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-0.3157};
        }
        return quadBondBondAdj;
    }

    protected double[] getEndBondTorPot(String atomOne, String atomTwo, String atomThree, String atomFour){
        double[] quadBondBondAdj = new double[1];
        if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C3A")){
            quadBondBondAdj = new double[]{-0.1185, 6.3204, 0.0000,-0.1185, 6.3204, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{ 0.0000,-0.6918, 0.0000, 0.0000, 0.2421, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.0000,-6.8958, 0.0000, 0.0000,-0.4669, 0.0000};
        }else if (atomOne.equals("C4") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.0000,-1.7970, 0.0000, 0.0000,-0.4879, 0.0000};
        }else if (atomOne.equals("H1") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.0000,-0.6890, 0.0000, 0.0000,-0.6890, 0.0000};
        }else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-0.5835, 1.1220, 0.3978, 1.3997, 0.7756, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{-0.0732, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.2486, 0.2422,-0.0925, 0.0814, 0.0591, 0.2219};
        }else if (atomOne.equals("H1") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.2130, 0.3120, 0.0777, 0.2130, 0.3120, 0.0777};
        }
        return quadBondBondAdj;
    }

    protected double[] getMidBondTorPot(String atomOne, String atomTwo, String atomThree, String atomFour){
        double[] quadBondBondAdj = new double[1];
        if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C3A")){
            quadBondBondAdj = new double[]{ 27.5989,-2.3120, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{ 0.0000, 9.1792, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{0.0000,-1.1521, 0.0000};
        }else if (atomOne.equals("C4") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.0000, 3.9421, 0.0000};
        }else if (atomOne.equals("H1") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.0000, 4.8228, 0.0000};
        }else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{ 0.0000, 0.0000, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-5.5679, 1.4083, 0.3010};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.0000, 0.0000, 0.0000};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{-17.7870,-7.1877, 0.0000};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-14.8790,-3.6581,-0.3138};
        }else if (atomOne.equals("H1") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-14.2610,-0.5322,-0.4864};
        }
        return quadBondBondAdj;
    }
    protected double[] getAngTorPot(String atomOne, String atomTwo, String atomThree, String atomFour){
        double[] quadBondBondAdj = new double[1];
        if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C3A")){
            quadBondBondAdj = new double[]{ 1.9767, 1.0239, 0.0000, 1.9767, 1.0239, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{ 0.0000, 3.8987, 0.0000, 0.0000,-4.4683, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.0000, 2.5014, 0.0000, 0.0000, 2.7147, 0.0000};
        }else if (atomOne.equals("C4") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.0000,-0.1242, 0.0000, 0.0000, 3.4601, 0.0000};
        }else if (atomOne.equals("H1") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.0000, 2.4501, 0.0000, 0.0000, 2.4501, 0.0000};
        }else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{ 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.2251, 0.6548, 0.1237, 4.6266, 0.1632, 0.0461};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{ 0.3886,-0.3139, 0.1389, 0.3886,-0.3139, 0.1389};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-0.2454, 0.0000,-0.1136, 0.3113, 0.4516,-0.1988};
        }else if (atomOne.equals("H1") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-0.8085, 0.5569,-0.2466,-0.8085, 0.5569,-0.2466};
        }
        return quadBondBondAdj;
    }

    protected double[] getAngAngTorPot(String atomOne, String atomTwo, String atomThree, String atomFour){
        double[] quadBondBondAdj = new double[1];
        if(atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("C4")){
            quadBondBondAdj = new double[]{ -14.4097};
        } else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-4.8141};
        } else if (atomOne.equals("C4") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 4.4444};
        }else if (atomOne.equals("H1") && atomTwo.equals("C3A") && atomThree.equals("C3A") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{ 0.3598};
        }else if (atomOne.equals("C3A") && atomTwo.equals("C3A") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-5.8888};
        }else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("C4")) {
            quadBondBondAdj = new double[]{-22.0450};
        } else if (atomOne.equals("C4") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-16.1640};
        } else if (atomOne.equals("H1") && atomTwo.equals("C4") && atomThree.equals("C4") && atomFour.equals("H1")) {
            quadBondBondAdj = new double[]{-12.5640};
        }
        return quadBondBondAdj;
    }
    public String[] getTwoAtomOrder (String atomNameOne, String atomNameTwo){
        String[] orderTwo = new String[2];
        int prefOne = getAtomTypePref(atomNameOne);
        int preftwo = getAtomTypePref(atomNameTwo);
        if(prefOne<preftwo || prefOne == preftwo){
            orderTwo[0] = atomNameOne;
            orderTwo[1] = atomNameTwo;
        } else {
            orderTwo[0] = atomNameTwo;
            orderTwo[1] = atomNameOne;
        }
        return orderTwo;
    }
    public String[] getThreeAtomOrder (String atomNameOne, String atomNameTwo, String atomNameThree){
        String[] orderTwo = new String[3];
        int prefOne = getAtomTypePref(atomNameOne);
        int prefthree = getAtomTypePref(atomNameThree);
        if(prefOne<prefthree || prefOne == prefthree){
            orderTwo[0] = atomNameOne;
            orderTwo[1] = atomNameTwo;
            orderTwo[2] = atomNameThree;
        } else {
            orderTwo[0] = atomNameThree;
            orderTwo[1] = atomNameTwo;
            orderTwo[2] = atomNameOne;
        }
        return orderTwo;
    }
    public String[] getFourAtomOrder (String atomNameOne, String atomNameTwo, String atomNameThree, String atomNameFour){
        String[] orderTwo = new String[4];
        int prefOne = getAtomTypePref(atomNameOne);
        int prefTwo = getAtomTypePref(atomNameTwo);
        int prefThree = getAtomTypePref(atomNameThree);
        int prefFour = getAtomTypePref(atomNameFour);
        if(prefTwo<prefThree || prefTwo ==prefThree){
            if(prefOne<prefFour || prefOne == prefFour){
                orderTwo[0] = atomNameOne;
                orderTwo[1] = atomNameTwo;
                orderTwo[2] = atomNameThree;
                orderTwo[3] = atomNameFour;
            }else {
                orderTwo[0] = atomNameFour;
                orderTwo[1] = atomNameThree;
                orderTwo[2] = atomNameTwo;
                orderTwo[3] = atomNameOne;
            }
        }else {
            orderTwo[0] = atomNameFour;
            orderTwo[1] = atomNameThree;
            orderTwo[2] = atomNameTwo;
            orderTwo[3] = atomNameOne;

        }
        return orderTwo;
    }
    public int getAtomTypePref(String elementName){
        int numPref =0;
        switch (elementName){
            case "C3A":
                numPref = 1;
                break;
            case "C4":
                numPref = 2;
                break;
            case "H1":
                numPref = 3;
                break;
        }
        return numPref;
    }

    public static void main(String[] args) {
        String atomOne = "H1";
        String atomTwo = "C3A";
        String atomThree ="C4";
        String atomFour ="C4";
       PDBReaderCOMPASS pdbReaderCOMPASS = new PDBReaderCOMPASS();
        String confName = "F://Avagadro//molecule//Ar";
        Vector centreMOP = new Vector3D(0.0,0.0, 0.0);
        pdbReaderCOMPASS.getSpecies(confName, false, false,centreMOP  );
       /*String[] orderTwo = pdbReaderCOMPASS.getTwoAtomOrder(atomOne, atomTwo);
        String[] orderTwoNew = pdbReaderCOMPASS.getTwoAtomOrder(atomTwo, atomThree);
        String[] orderTwoNewNew = pdbReaderCOMPASS.getTwoAtomOrder(atomThree, atomFour);
        //String[] orderThree = pdbReaderCOMPASS.getThreeAtomOrder(atomOne, atomTwo, atomThree);
       // String[] orderFour = pdbReaderCOMPASS.getFourAtomOrder(atomOne, atomTwo, atomThree, atomFour);
       // System.out.println(Arrays.toString(orderTwo) + Arrays.toString(orderTwoNew));
        //System.out.println(Arrays.toString(orderThree));
        double[] valueBondOne = pdbReaderCOMPASS.getQuadBondPot(orderTwo[0], orderTwo[1]);
        double[] valueBondTwo = pdbReaderCOMPASS.getQuadBondPot(orderTwoNew[0], orderTwoNew[1]);
        double[] valueBondThree = pdbReaderCOMPASS.getQuadBondPot(orderTwoNewNew[0], orderTwoNewNew[1]);
        System.out.println(Arrays.toString(orderTwo) + " " + Arrays.toString(orderTwoNew) + " " + Arrays.toString(orderTwoNewNew));
        System.out.println(Arrays.toString(valueBondOne) + " " + Arrays.toString(valueBondTwo) + " " + Arrays.toString(valueBondThree));
    }

}*/
