package etomica.potential.UFF;


import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.chem.elements.*;
import etomica.molecule.IMolecule;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class PDBReaderReplicaNew {
    int a, b, c;
    String atomA, atomB, atomC, atomTypeA, atomTypeB, atomTypeC;
    public Map<Integer,Vector> positions = new HashMap<>();
    public Map<Integer, String> atomMap = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectivity = new ArrayList<>();
    public HashMap<Integer, String> atomMapModified = new HashMap<>();
    public Map<Integer, Integer> coordinationMapOutput  = new HashMap<>();
    public Map<String, AtomType> typeMap = new HashMap<>();
    public Map<String, AtomType> typeMapNew = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectivityModified = new ArrayList<>();
    Map<Integer, String> modifiedAtomIdentifierMap = new HashMap<>();
    public List<int[]> listOfInversions = new ArrayList<>();
    public ArrayList<Integer> bondList = new ArrayList<>();
    public Map<Integer, AtomType> atomIdentifierMap = new HashMap<>();
    public ArrayList<Integer> bondsNum = new ArrayList<>();
    public List<int[]> intListArray = new ArrayList<>();
    public List<int[]> tripletsList = new ArrayList<>();
    public List<int[]> quadruplets = new ArrayList<>();
    public List<int[]> dupletsSorted = new ArrayList<>();
    public List<int[]> tripletsSorted = new ArrayList<>();
    public List<int[]> quadrupletsSorted = new ArrayList<>();
    public Map<String, Element> elementReaderMap = new HashMap<>();
    public Map<String, AtomType> elementReceiverMap = new HashMap<>();
    public Map<AtomType, AtomType> elementReceiverMapNew= new HashMap<>();
    ISpecies species;
    int connect = 0;

    public void readPDBFile(String confName) {
        String fileName = confName+".pdb";
        FileReader fileReader;
        Map<String, AtomType> typeMap = new HashMap<>();
        ArrayList<Integer> currentAtomList = null;
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            String line ;

            while ((line = bufReader.readLine()) != null) {

                parseLineReader(line, typeMap, atomMap, positions);

                if (line.startsWith("CONECT")) {
                    String[] parts = line.trim().split("\\s+");
                    int atomNumber = Integer.parseInt(parts[1]);
                    if (currentAtomList == null || atomNumber != currentAtomList.get(0)) {
                        // start a new list for the current atom
                        currentAtomList = new ArrayList<>();
                        connectivity.add(currentAtomList);
                        currentAtomList.add(atomNumber);
                        //System.out.println(currentAtomList + "up");
                    }
                    for (int i = 2; i < parts.length; i++) {
                        currentAtomList.add(Integer.parseInt(parts[i]));
                        //System.out.println(currentAtomList + "down");
                    }
                    connect++;
                }
            }
            /*if(connect == 0){
                System.out.println("Nobel Gas");
                ArrayList<Integer> nobel = new ArrayList<>();
                nobel.add(0, 0);
                connectivity.add(nobel);
            }*/
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }

        //String  molecularWeightString = String.format("%.4f", molecularWeight);
        // System.out.println( "Molecular Weight: "+molecularWeightString);
        System.out.println(typeMap + "typeMap");
        System.out.println(connectivity + " connectivity");
        System.out.println(atomMap+" atomMap");
        System.out.println("Building done");
        //COM(positions);
        System.out.println("After");
        System.out.println(positions);
        //System.exit(1);
    }

    protected void parseLineReader(String line, Map<String, AtomType > typeMap, Map<Integer, String> atomMap, Map<Integer, Vector> positions) {
        line = line.trim();
        if (line.length() < 6) {
            return;
        }
        if (line.substring(0, 6).equals("HETATM") || line.substring(0,4).equals("ATOM")) {
            //coordinates of atom and create atom
            double x = Double.parseDouble(line.substring(30, 38));
            double y = Double.parseDouble(line.substring(38, 46));
            double z = Double.parseDouble(line.substring(46, 54));
            Vector positn = Vector.of(x, y, z);
            String symbol = line.substring(12, 16).trim();
            int atomNumber = Integer.parseInt(line.substring(8,11).trim());
            System.out.println(atomNumber + " " + positn);
            positions.put(atomNumber, positn);
           /* if (typeMap.containsKey(symbol)) {
                //type = typeMap.get(symbol);
                typeElement = putElements(line.substring(76, 78).trim());
                type = AtomType.element(typeElement);
               // getCOM(typeElement,positn);
                //System.out.println(type + "contains");
            } else {
                typeElement = putElements(line.substring(76, 78).trim());
                type = AtomType.element(typeElement);
                //getCOM(typeElement,positn);
                //System.out.println(type + " newtype");
                typeMap.put(symbol, type);
            }*/
            if (atomMap.containsKey(atomNumber)){
                System.out.println("Error");
            } else {
                // tempMap.put(atomNumber, typeElement);
                atomMap.put(atomNumber, symbol);
            }
        }
    }

    public Map<Integer, String> getAtomMap(){
        return atomMap;
    }
    public Map<Integer,String>getAtomModifiedMap(){
        return atomMapModified;
    }

    public ISpecies getSpecies (String confName){
        double massSum = 0;
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        SpeciesBuilder speciesBuilderNewMod =  new SpeciesBuilder(Space3D.getInstance());
        AtomType typeNew;
        readPDBFile(confName);
        System.out.println(connectivity);
        System.out.println(atomMap);
        if(connect !=0){
            ArrayList<Integer> bondList =setBondList(connectivity, atomMap);
           System.out.println(bondList + " bondList");
        }
        System.out.println(atomMap+ "  " + atomMap.size());
        Map<Integer, String> atomIdentifierMapMod = atomIdentifierMapModified(connectivity, atomMap);
        // System.out.println(atomIdentifierMapModified + "atomIdentifierMapModified");
        ArrayList<ArrayList<Integer>> connectedAtoms =getConnectivity();
        //System.out.println(connectedAtoms+ ": connectedAtom");
        ArrayList<ArrayList<Integer>> connectivityModified = getconnectivityModified(connectedAtoms);
        // System.out.println(connectivityModified+ ": connectedAtomModified" );
        Map<Integer,String> atomMap = getAtomMap(connectedAtoms);
        // System.out.println(atomMap + ": atomMap");
        HashMap<Integer, String> atomMapModified =getatomMapModified(atomMap);
        // System.out.println(atomMapModified + ": atomMapModified");

        System.out.println(connectedAtoms +" O");
        System.out.println(connectivityModified+" O");
        System.out.println(atomMapModified+" O");
        System.out.println(atomIdentifierMapMod+" O");
        if(connect !=0){
            bondsNum = bondsAmongAtoms();
            System.out.println(bondsNum + " bond Amongatoms");
            List<int[]> duplets = getBondedAtomList(connectivityModified);
            List<int[]> listOfAngleModified = getAngleList(connectivityModified);
            dupletsSorted= bondSorter(duplets, atomIdentifierMapMod);
            //System.out.println(Arrays.deepToString(dupletsSorted.toArray()) + " duplets Sorted");
            //System.out.println(Arrays.deepToString(duplets.toArray())+ ": listOfBonds");
            tripletsSorted=angleSorter(listOfAngleModified, atomIdentifierMapMod);
            List<int[]> listOfTorsionModified = getTorsionList(connectivity);
            //System.out.println( Arrays.deepToString(listOfTorsionModified.toArray()) + "torsionModified");
            quadrupletsSorted=torsionSorter(listOfTorsionModified, atomIdentifierMapMod);
            // System.out.println(Arrays.deepToString(quadrupletsSorted.toArray()) + "quadrupletsSorted");
            Map<String[],List<int[]>> torsionTypesMap= idenTorsionTypes(quadrupletsSorted, atomIdentifierMapMod);
            aromaticIdentifier(torsionTypesMap);
            aromaticOtherElementFormer(connectivityModified, modifiedAtomIdentifierMap);
            listOfInversions = idenInversions(tripletsSorted);
            //System.out.println(Arrays.deepToString(listOfInversions.toArray()) + " in Main" );
        }
        System.out.println(positions + " Here positions");
        Space space = Space3D.getInstance();
        Vector center = space.makeVector();
        Vector dr = Vector.d(center.getD());
        System.out.println(atomIdentifierMap);
        for(int i = 0; i < atomIdentifierMap.size(); i++) {
            String symbol = String.valueOf(atomIdentifierMap.get(i));
            int startIndex = symbol.indexOf("[") + 1;
            int endIndex = symbol.indexOf("]");
            String nameNew = symbol.substring(startIndex, endIndex);
            // System.out.println(atomIdentifierMap.get(i));
            AtomType newName = returnElement(nameNew);
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
        System.out.println(typeMapNew + " typeMapNew");
        species= speciesBuilderNew.build();
        System.out.println(species.getMass() + " first");
        IMolecule molecule = species.makeMolecule();
        IAtomList children = molecule.getChildList();
        int nAtoms = children.size();
        for (int i = 0; i < nAtoms; i++) {
            IAtom a = children.get(i);
            // System.out.println(a.getPosition() + " "+ i);
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
        // System.out.println(center + " 1 ");
        // System.out.println("Part 1");
        for(int i=0; i<atomIdentifierMap.size(); i++){
            IAtom a = children.get(i);
            String symbol = String.valueOf(atomIdentifierMap.get(i));
            int startIndex = symbol.indexOf("[") + 1;
            int endIndex = symbol.indexOf("]");
            String nameNew = symbol.substring(startIndex, endIndex);
            typeNew = typeMapNew.get(nameNew);
            Vector v = a.getPosition();
            v.ME(center);
            // System.out.println(v);
            speciesBuilderNewMod.addAtom(typeNew,v, "" );
        }
        species = speciesBuilderNewMod.build();

        //System.out.println(species.getMass() + " Second");
        return species;
    }
    public List<int[]> getDupletesSorted(){
        return dupletsSorted;
    }
    public List<int[]> getAnglesSorted(){
        return tripletsSorted;
    }
    public List<int[]> getTorsionSorted(){
        return quadrupletsSorted;
    }


    public void clearEverything(){
        positions.clear();
        atomMap.clear();
        connectivity.clear();
        atomMapModified.clear();
        coordinationMapOutput.clear();
        typeMapNew.clear();
        connectivityModified.clear();
        modifiedAtomIdentifierMap.clear();
        listOfInversions.clear();
        bondList.clear();
        atomIdentifierMap.clear();
        bondsNum.clear();
        intListArray.clear();
        tripletsList.clear();
        quadruplets.clear();
        dupletsSorted.clear();
        tripletsSorted.clear();
        quadrupletsSorted.clear();
    }

    public void aromaticIdentifier( Map<String[],List<int[]>> torsionTypesMap){
        //System.out.println(atomIdentifierMap);
        //System.out.println(connectivityModified);
        for (Map.Entry<String[], List<int[]>> entry : torsionTypesMap.entrySet()) {
            String[] torsionType = entry.getKey();
            List<int[]> torsions = entry.getValue();
            //System.out.println(Arrays.toString(torsionType) +" torsion");
            if(torsionType[0].equals("C_2") && torsionType[3].equals("C_2") && torsionType[1].equals("C_2") && torsionType[2].equals("C_2")){
                // System.out.println("entered loop " + Arrays.deepToString(torsions.toArray()) );
                for(int i=0; i<torsions.size();i++){
                    int[] torsionIndividual = torsions.get(i);
                    System.out.println(Arrays.toString(torsionIndividual));
                    int torsionIndividualOne = torsionIndividual[1];
                    int torsionIndividualTwo = torsionIndividual[2];
                    //  String torsionIndividualElementOne = modifiedAtomIdentifierMap.get(torsionIndividualOne);
                    // String torsionIndividualElementTwo = modifiedAtomIdentifierMap.get(torsionIndividualTwo);
                    // if(torsionIndividualElementOne.equals("C_Ar") && torsionIndividualElementTwo.equals("C_Ar")){
                    //    i++;
                    //} else {
                    int[] hexAromaticNew = new int[6];
                    System.out.println(torsionIndividualOne + " " +torsionIndividualTwo);
                    ArrayList<Integer> connectivityFirst = connectivityModified.get(torsionIndividual[0]);
                    ArrayList<Integer> connectivityLast = connectivityModified.get(torsionIndividual[3]);
                    int start = -1;
                    for(int j=1; j<connectivityFirst.size(); j++){
                        // System.out.println(connectivityFirst +" connectivity First");
                        int atomNum = connectivityFirst.get(j);
                        String atomName = modifiedAtomIdentifierMap.get(atomNum);
                        // System.out.println(atomNum + " atomName "+ atomName);
                        if( (atomName.equals("C_2") || atomName.equals("C_Ar")) && atomNum != torsionIndividualOne){
                            start = atomNum;
                        }
                        System.out.println(start);
                        // System.out.println(modifiedAtomIdentifierMap.get(connectivityFirst.get(j)) + " " +connectivityFirst.get(j));
                    }
                    int end = -1;
                    for(int j=1; j<connectivityLast.size(); j++) {
                        // System.out.println(connectivityLast +" connectivityLast");
                        int atomNum = connectivityLast.get(j);
                        String atomName = modifiedAtomIdentifierMap.get(atomNum);
                        if( (atomName.equals("C_2") || atomName.equals("C_Ar")) && atomNum != torsionIndividualTwo) {
                            end = atomNum;
                        }
                        System.out.println(end);
                        // System.out.println(modifiedAtomIdentifierMap.get(connectivityFirst.get(j)) + " " +connectivityFirst.get(j));
                    }
                    System.out.println(start + " start");
                    if(start > -1 && end > -1){
                        ArrayList<Integer> startOneConnect = connectivityModified.get(start);
                        // ArrayList<Integer> endOneConnect = connectivityModified.get(endOne);
                        ArrayList<Integer> aromaticRing = new ArrayList<>();
                        for(int k =0; k<startOneConnect.size(); k++){
                            if(end == startOneConnect.get(k)){
                                aromaticRing.add(start);
                                aromaticRing.add(torsionIndividual[0]);
                                aromaticRing.add(torsionIndividual[1]);
                                aromaticRing.add(torsionIndividual[2]);
                                aromaticRing.add(torsionIndividual[3]);
                                aromaticRing.add(end);
                                // System.out.println(aromaticRing + " ring in loop");
                                aromaticConverter( aromaticRing, atomIdentifierMap );
                            }
                        }
                        System.out.println( Arrays.toString(torsionIndividual) + " " + connectivityLast + " "  );
                    }

                }
            } else {
                //  System.out.println("Non Aromatic");
            }
        }
    }

    public void aromaticConverter(ArrayList<Integer> ring, Map<Integer, AtomType> atomIdentifierMap){
        // System.out.println(atomIdentifierMap +" atomIdentifierMap");
        //System.out.println(atomMap + " atomMap");
        // System.out.println(atomMapModified + " atommapModified");
        // System.out.println(modifiedAtomIdentifierMap + " modifiedAtomIdentifier");
        // System.out.println(ring +" ring");
        for(int i =0; i< ring.size(); i++){
            IElement Carbon_Ar = null;
            AtomType C_Ar = new AtomType(Carbon_Ar, "C_Ar");
            System.out.println(ring.get(i));
            atomIdentifierMap.put(ring.get(i), C_Ar);
            modifiedAtomIdentifierMap.put(ring.get(i), "C_Ar");
            System.out.println(i +" " + atomIdentifierMap.get(i) );
        }
        // System.out.println(modifiedAtomIdentifierMap + " map modified");
        // System.out.println(atomIdentifierMap + " atomIdentifierMap");
    }
    private void aromaticOtherElementFormer(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> modifiedAtomIdentifierMap) {
        for(int i=0; i<connectivityModified.size(); i++){
            ArrayList<Integer> connect = connectivityModified.get(i);
            int atomNumZero = connectivityModified.get(i).get(0);
            String atomName = modifiedAtomIdentifierMap.get(atomNumZero);
            if(atomName.equals("C_Ar")){
                for(int j =1; j<connect.size(); j++ ){

                }
                /* if (atomName.equals("O_3") || atomName.equals("O_2")) {
                IElement Oxygen_Ar = null;
                AtomType O_Ar = new AtomType(Oxygen_Ar, "O_Ar");
                //System.out.println("Oxy");
                modifiedAtomIdentifierMap.put(i, "O_Ar");
                atomIdentifierMap.put(i, O_Ar);
            } else if (atomName.equals("N_3") || atomName.equals("N_2")) {
                IElement Nitrogen_Ar = null;
                AtomType N_Ar = new AtomType(Nitrogen_Ar, "N_Ar");
                //System.out.println("Nitro");
                modifiedAtomIdentifierMap.put(i, "N_Ar");
                atomIdentifierMap.put(i, N_Ar);
            }*/
            } else {
                break;
            }


        }
        // System.out.println(atomIdentifierMap +" atomIdentifierMap Final");
        // System.out.println(modifiedAtomIdentifierMap +" modified AtomMap final");
    }


    public ArrayList<Integer> bondsAmongAtoms(){
        ArrayList<Integer> bondListSatisfied = new ArrayList<>();
        bondListSatisfied.addAll(bondList);
        //System.out.println(bondListSatisfied + " bondlist satisfied");
        //System.out.println(bondList + " bondlist");
        ArrayList<int[]> bondsAmongAtoms = new ArrayList<>();

        // System.out.println("bonds Among Atoms");

        System.out.println(connectivityModified);
        System.out.println(atomIdentifierMap);
        for(int i =0; i<connectivityModified.size(); i++){
            int num = connectivityModified.get(i).get(0);
            System.out.println(connectivityModified.get(i) +" "+ atomIdentifierMap.get(i));
            for(int j=1; j<connectivityModified.get(i).size(); j++){
                int val = connectivityModified.get(i).get(j);
                if(num<val){
                    //System.out.println(atomIdentifierMap.get(i) + " " + bondListSatisfied.get(i) + " " +val + " " +  bondListSatisfied.get(val) );
                    if(bondListSatisfied.get(i) == 2 && bondListSatisfied.get(val) ==2){
                        //System.out.println("2");
                        bondListSatisfied.set(i, 1);
                        bondListSatisfied.set(val,1);
                        bondsNum.add(2);
                    } else {
                        System.out.println("1");
                        bondsNum.add(1);
                    }
                }
                //System.out.println(bondsNum);
            }
        }
        return bondsNum;
    }
    public ArrayList<Integer> getBonds(){
        return bondsNum;
    }

    public List<int[]> getListofInversions(){
        return listOfInversions;
    }


    public List<int[]> idenInversions(List<int[]> angleList){
        //List<int[]> inversionList = new ArrayList<>();
        System.out.println(atomIdentifierMap);
        //  System.out.println(angleList.size());
        int k = 0;
        for(int i=0; i<angleList.size(); i++){
            int[] arr= new int[4];
            ArrayList<Integer> arrayList = new ArrayList<>();
            int[] angle = angleList.get(i);
            a = angle[0];
            b = angle[1];
            c = angle[2];
            atomA = String.valueOf(atomIdentifierMap.get(a));
            atomB = String.valueOf(atomIdentifierMap.get(b));
            atomC = String.valueOf(atomIdentifierMap.get(a));
            atomTypeA = atomA.substring(9, atomA.length() - 1);
            atomTypeB = atomB.substring(9, atomB.length() - 1);
            atomTypeC = atomC.substring(9, atomC.length() - 1);
            if(atomTypeA.equals("C_Ar") && atomTypeB.equals("C_Ar") && atomTypeC.equals("C_Ar")){
                arrayList = connectivityModified.get(b);
                arr =  arrayList.stream().mapToInt(j -> j).toArray();
                if(k%3 == 0){
                    System.out.println(Arrays.toString(arr));
                    System.out.println(Arrays.deepToString(listOfInversions.toArray()));
                    listOfInversions.add(arr);
                }
                k++;
            }
            //  System.out.println(Arrays.toString(angle) +" " + a +" "+b +" "+c);
        }
        //System.out.println();
        return listOfInversions;
    }
    public void makeBondArray(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> modifiedAtomIdentifierMap, ArrayList<Integer> bondList){
        System.out.println("MakeBondArray");
        System.out.println(connectivityModified);
        System.out.println(modifiedAtomIdentifierMap);
        int[][] bondArray = new int[modifiedAtomIdentifierMap.size()][modifiedAtomIdentifierMap.size()];
        for(int i =0; i<modifiedAtomIdentifierMap.size(); i++){
            ArrayList<Integer> connectivity_array = connectivityModified.get(i);
            System.out.println(connectivity_array);
            int firstElement = connectivity_array.get(0);
            int bondListA = bondList.get(firstElement);
            int enter = 0;
            for(int j=1; j<connectivity_array.size(); j++){
                int particularElement = connectivity_array.get(j);
                //System.out.println(connectivity_array.get(j) + " " + firstElement);
                int bondListB = bondList.get(particularElement);
                if(bondListA == 2 && bondListB == 2 && enter < 1){
                    bondArray[firstElement][particularElement] = 2;
                    enter ++;
                } else {
                    bondArray[firstElement][particularElement] = 1;
                }

            }
        }
        System.out.println(Arrays.deepToString(bondArray) + " bondArray");
    }

    public Map<Integer,AtomType> getAtomIdentifierMapModified(){
        return atomIdentifierMap;
    }

    public Map<Integer,String> atomIdentifierMapModified (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified) {
        Map<Integer, AtomType> atomIdentifierMap = atomIdentifier(connectivityModified, atomMapModified);

        for (Map.Entry<Integer, AtomType> entry : atomIdentifierMap.entrySet()) {
            String value = entry.getValue().toString();

            value = value.replace("AtomType[", "").replace("]", "");
            modifiedAtomIdentifierMap.put(entry.getKey(), value);
        }
        System.out.println(modifiedAtomIdentifierMap + " modified Identifier");
        return modifiedAtomIdentifierMap;
    }
    public ArrayList<Integer> getBondList (ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> atomMap){
        return bondList;
    }
    public Map<Integer,String> getatomIdentifierMapModified(){
        return modifiedAtomIdentifierMap;
    }
    public ArrayList<Integer> setBondList(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        ArrayList<Integer> bondList = singleDoubleBondIdentifier(connectivityModified, atomMapModified);
        return bondList;
    }
    public ArrayList<Integer> singleDoubleBondIdentifier(ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> atomMap) {
        //System.out.println("Already visited");
        int valenceCarbon = 4, valenceOxygen = 2, valenceHydrogen = 1, valenceNitrogen =0, valenceSulfur = 0, valenceHalide =1;
        int bondRequired = 0;
        ArrayList<Integer> nitrogenBonds;
        ArrayList<Integer> sulfurBonds;
        System.out.println(atomMap);
        for (int i = 0; i < connectivity.size(); i++) {
            int atomArraySize = connectivity.get(i).size();
            // System.out.println(atomArraySize + " start atomarraysize" +  connectivity.get(i));
            String atomName = atomMap.get(connectivity.get(i).get(0));
            if (atomName.equals("H")) {
                if (valenceHydrogen == atomArraySize - 1) {
                    //System.out.println("The hydrogen atom " + i + " is satisfied with one bond");
                    bondList.add(1);
                } else {
                    //bondRequired = valenceHydrogen - (atomArraySize - 1);
                    // System.out.println("The hydrogen atom " + i + " is UNSATISFIED with " + bondRequired + " bond required");
                    bondList.add(0);
                }
            }
            if (atomName.equals("CL") ||atomName.equals("BR") || atomName.equals("F")) {
                if (valenceHydrogen == atomArraySize - 1) {
                    // System.out.println("The hydrogen atom " + i + " is satisfied with one bond");
                    bondList.add(1);
                } else {
                    //bondRequired = valenceHydrogen - (atomArraySize - 1);
                    //System.out.println("The hydrogen atom " + i + " is UNSATISFIED with " + bondRequired + " bond required");
                    bondList.add(0);
                }
            }

            if (atomName.equals("O")) {
                if (valenceOxygen == atomArraySize - 1) {
                    //   System.out.println("The oxygen atom " + i + " is satisfied with one bond");
                    bondList.add(1);
                } else if (valenceOxygen == atomArraySize ) {
                    bondList.add(2);
                } else {
                    // bondRequired = valenceOxygen - (atomArraySize - 1);
                    // System.out.println("The oxygen atom " + i + "is UNSATISFIED with " + bondRequired + " bond required");
                    int numNitro=0;
                    for(int j =1; j<connectivity.get(i).size();j++){
                        String atomConnectedtoOxygen = atomMap.get(connectivity.get(i).get(j));
                        if(atomConnectedtoOxygen.equals("N")){
                            numNitro++;
                        }
                    }
                }
            }

            if (atomName.equals("C")) {
                //System.out.println(atomArraySize);
                if(atomArraySize == 3){
                    bondList.add(4);
                }
                if (valenceCarbon == atomArraySize - 1) {
                    // System.out.println("The carbon atom " + i + " is satisfied with one bond");
                    bondList.add(1);
                } else {
                    //bondRequired = valenceCarbon - (atomArraySize - 1);
                    //   System.out.println("The carbon atom " + i + " is UNSATISFIED with " + bondRequired + " bond required");
                    bondList.add(0);
                }
            }

            if (atomName.equals("N")){
                nitrogenBonds = connectivity.get(i);
                // System.out.println(nitrogenBonds);
                int numCarb=0, numNitro=0, numOxy=0, numHydro=0, numMetal=0;
                for(int j =1; j<connectivity.get(i).size();j++){
                    String atomConnectedtoNitrogen = atomMap.get(connectivity.get(i).get(j));
                    if (atomConnectedtoNitrogen.equals("C")){
                        numCarb++;
                    } else if(atomConnectedtoNitrogen.equals("H")){
                        numHydro++;
                    } else if(atomConnectedtoNitrogen.equals("O")){
                        numOxy++;
                    } else if (atomConnectedtoNitrogen.equals("N")) {
                        numNitro++;
                    } else {
                        numMetal++;
                    }
                }
                //System.out.println(numCarb + ":numCarb " + numHydro +":numHydro " +  numOxy +":numOxy " +  numNitro + ":numNitro " + numMetal+":numMetal" );
                if (numMetal >1){
                    // valenceNitrogen =4;
                    if( atomArraySize < 5){
                        bondList.add(0);
                    } else {
                        bondList.add(1);
                    }
                } else if (numOxy ==2) {
                    // valenceNitrogen=4;
                    if( atomArraySize < 5){
                        bondList.add(0);
                    } else {
                        bondList.add(1);
                    }
                } else {
                    // valenceNitrogen =3;
                    if(atomArraySize < 4){
                        bondList.add(0);
                    } else {
                        bondList.add(1);
                    }
                }
                //System.out.println(valenceNitrogen + "valenceNitrogen");
            }
            if(atomName.equals("S")){
                sulfurBonds = connectivity.get(i);
                System.out.println(sulfurBonds);
                int numCarb=0, numOxy=0;
                for(int j =1; j<connectivity.get(i).size();j++){
                    String atomConnectedtoNitrogen = atomMap.get(connectivity.get(i).get(j));
                    if (atomConnectedtoNitrogen.equals("C")){
                        numCarb++;
                    } else if(atomConnectedtoNitrogen.equals("O")) {
                        numOxy++;
                    }
                }
                if ((numCarb ==2 && numOxy ==0)||(numCarb ==1 && numOxy != 3)){
                    valenceSulfur = 2;
                } else if ((numCarb ==2 && numOxy ==2) ||(numCarb ==1 && numOxy ==3)){
                    valenceSulfur = 4;
                }
                // System.out.println(valenceSulfur + "valencesulfur");
            }

        }
        // System.out.println("Bonds are satisfied or not: "+ bondList);
        System.out.println(atomMap);
        ArrayList<Integer> connectivityArraylist = new ArrayList<>();
        ArrayList<Integer> connectivityArrayElementBondList = new ArrayList<>();
        int index;
        int e=0;
        while (e<=5){
            // System.out.println("Reached after while loop");
            for(int i =0; i<bondList.size(); i++){
                //  System.out.println(bondList.get(i) + " bondlist element " + i   );
                if (bondList.get(i) ==0){
                    // System.out.println((i+1) + " " + bondList.get(i) + "i ");
                    //System.out.println(connectivity.get(i) + " i" + i + "connectivity here");
                    connectivityArraylist = connectivity.get(i);
                    //System.out.println(connectivityArraylist + " connectivityarraylist");
                    for(int j =0; j<connectivityArraylist.size(); j++){
                        connectivityArrayElementBondList =connectivity.get(connectivityArraylist.get(j)-1);
                        index = connectivity.indexOf(connectivityArrayElementBondList);
                        if(bondList.get(index) == 0) {
                            // System.out.println(connectivityArrayElementBondList + " connectivityArrayElementBondList " + " " + bondList.get(index) + " bondlist element " + j + 1);
                            // System.out.println(connectivityArraylist + "Arraylist");
                            // System.out.println(connectivityArrayElementBondList + " arrayelmentlist");
                            if(connectivityArraylist != connectivityArrayElementBondList){
                                bondList.set(i,1);
                                bondList.set(index, 1);
                                //System.out.println("The atoms linked are" + connectivityArraylist+ " " +(i+1) + "  is connected to " + connectivityArrayElementBondList +" index: "+ (index+1) + " are double bonded");
                                bondList.set(i, 2);
                                bondList.set(index, 2);
                                j++;
                            }
                        }

                    }
                    //System.out.println(bondList + "bondlist");
                }
            }
            e++;
        }
        if(bondList.size() < atomMap.size()){
            bondList.add(0);
        }
        return bondList;
    }
    public int priorityMapGenerator(String atomType){
        Map<String,Integer> priorityMap = new HashMap<>();
        priorityMap.put("Zr",1);
        priorityMap.put("C_3",2);
        priorityMap.put("C_3p",2);
        priorityMap.put("C_Arp",2);
        priorityMap.put("C_2p",3);
        priorityMap.put("C_1p",4);
        priorityMap.put("H",5);
        priorityMap.put("H_p",5);
        priorityMap.put("Ph",6);
        priorityMap.put("N_3p",7);
        priorityMap.put("N_2p",8);
        priorityMap.put("N_1p",9);
        priorityMap.put("O_3p",10);
        priorityMap.put("O_2p",11);
        priorityMap.put("O_1p",12);
        priorityMap.put("I",13);
        priorityMap.put("Br",14);
        priorityMap.put("S_3p", 6);
        priorityMap.put("S_2p", 6);
        priorityMap.put("S_1p", 6);
        priorityMap.put("Cl",15);
        priorityMap.put("F",16);
        priorityMap.put("Co", 16);
        return priorityMap.get(atomType);
    }

    public AtomType returnElement(String elementName){
        elementReceiverMap.put("C_3p", new AtomType(Carbon.INSTANCE, "C_3p"));
        elementReceiverMap.put("C_2p",  new AtomType(Carbon.INSTANCE, "C_2p"));
        elementReceiverMap.put("C_Arp",  new AtomType(Carbon.INSTANCE, "C_Arp"));
        elementReceiverMap.put("C_1p",  new AtomType(Carbon.INSTANCE, "C_1p"));
        elementReceiverMap.put("H_p",  new AtomType(Hydrogen.INSTANCE, "H_p"));
        elementReceiverMap.put("O_3p", new AtomType(Oxygen.INSTANCE, "O_3p"));
        elementReceiverMap.put("O_2p", new AtomType(Oxygen.INSTANCE, "O_2p"));
        elementReceiverMap.put("N_3p", new AtomType(Nitrogen.INSTANCE, "N_3p"));
        elementReceiverMap.put("N_2p", new AtomType(Nitrogen.INSTANCE, "N_2p"));
        elementReceiverMap.put("N_1p", new AtomType(Nitrogen.INSTANCE, "N_1p"));
        elementReceiverMap.put("S_3p", new AtomType(Sulfur.INSTANCE, "S_3p"));
        elementReceiverMap.put("S_2p",  new AtomType(Sulfur.INSTANCE, "S_2p"));
        elementReceiverMap.put("Ar",  new AtomType(Argon.INSTANCE, "Ar"));
        elementReceiverMap.put("He",  new AtomType(Helium.INSTANCE, "He"));
        elementReceiverMap.put("Ne",  new AtomType(Neon.INSTANCE, "Ne"));
        elementReceiverMap.put("Co",  new AtomType(Cobalt.INSTANCE, "Co"));
        elementReceiverMap.put("Cl",  new AtomType(Chlorine.INSTANCE, "Cl"));
        return elementReceiverMap.get(elementName);
    }

    protected Map<Integer, AtomType> atomIdentifier(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        //System.out.println(atomMapModified + " atommapModified here");
        //int counter =0;
        // ArrayList <Integer> atomNumbers = new ArrayList<>();
        if(connectivity.size() == 0){
            String element = atomMap.get(1);
            if(element.equals("AR")){
                AtomType Ar = new AtomType(Argon.INSTANCE, "Ar");
                atomIdentifierMap.put(0, Ar);
            } else if (element.equals("NE")) {
                AtomType Ne = new AtomType(Neon.INSTANCE, "Ne");
                atomIdentifierMap.put(0, Ne);
            } else {
                AtomType He = new AtomType(Helium.INSTANCE, "He");
                atomIdentifierMap.put(0, He);
            }
        }
        //System.out.println(connectivityModified + "connectivityModified");
        for (int i = 0; i < connectivityModified.size(); i++) {
            Integer retriveArrayFirstElementNumber = connectivityModified.get(i).get(0);
            //System.out.println(atomMapModified + "atomMap Modified");
            //System.out.println(connectivityModified + " connectivityModified");
            String retriveArrayFirstElementName = atomMapModified.get(retriveArrayFirstElementNumber);
            //System.out.println(retriveArrayFirstElementName + " :retriveFirstName");
            //System.out.println(atomMapModified + "atomMapModified");
            //Carbon atoms
            if (retriveArrayFirstElementName.equals("C")) {
                int arrayListSize = connectivityModified.get(i).size();
                //Alkane checker
                if (arrayListSize == 5) {
                    AtomType C_3p = new AtomType(Carbon.INSTANCE, "C_3p");
                    //System.out.println("The atom " + (i) +" is C_3" );
                    atomIdentifierMap.put(i, C_3p);
                }

                if (arrayListSize == 4) {
                    AtomType C_2p = new AtomType(Carbon.INSTANCE, "C_2p");
                    //System.out.println("The atom " + (i) +" is C_2 " );
                    atomIdentifierMap.put(i, C_2p);
                }
                if (arrayListSize == 3 || arrayListSize ==2) {
                    AtomType C_1p = new  AtomType(Carbon.INSTANCE, "C_1p");
                    //System.out.println("The atom " + (i)+" is C_1 "  );
                    atomIdentifierMap.put(i, C_1p);
                }
            } else if (retriveArrayFirstElementName.equals("CL")) {
                AtomType Cl = new AtomType(Chlorine.INSTANCE, "Cl");
                //System.out.println("The atom " + (i)+" is H " );
                atomIdentifierMap.put(i, Cl);
            } else if (retriveArrayFirstElementName.equals("H")) {
                AtomType H_p = new AtomType(Hydrogen.INSTANCE, "H_p");
                //System.out.println("The atom " + (i)+" is H " );
                atomIdentifierMap.put(i, H_p);
            }
            else if (retriveArrayFirstElementName.equals("O")) {
                int arrayListSize = connectivityModified.get(i).size();
                if(arrayListSize == 2){
                    AtomType O_2p = new AtomType(Oxygen.INSTANCE, "O_2p");
                    //System.out.println("The atom " + (i)+" is O_2 " );
                    atomIdentifierMap.put(i, O_2p);
                } else if (arrayListSize == 3) {

                    AtomType O_3p =new AtomType(Oxygen.INSTANCE, "O_3p");
                    //System.out.println("The atom " + (i ) + " is O_3 ");
                    atomIdentifierMap.put(i, O_3p);
                } else {
                    AtomType O_Ar = new AtomType(Oxygen.INSTANCE, "O_Ar");
                    //System.out.println("The atom " + (i ) + " is O_Ar ");
                    atomIdentifierMap.put(i, O_Ar);
                }
            } else if (retriveArrayFirstElementName.equals("N")){
                int arrayListSize = connectivityModified.get(i).size();
                int  oxygenCounter =0;
                //Types are:
                //-NH2, =N-, =N- -(+), -N=N-, O-(-)N(+)=O
                //System.out.println("It is Nitrogen");
                for (int j = 0; j< arrayListSize; j++){
                    String connectedElement = atomMapModified.get(connectivityModified.get(i).get(j));
                    if(connectedElement.equals("O")){
                        oxygenCounter++;
                    }
                }
                if(arrayListSize ==3){
                    AtomType N_1p = new AtomType(Nitrogen.INSTANCE, "N_1p");
                    //System.out.println("The atom " + (i) +" is N_1 " );
                    atomIdentifierMap.put(i, N_1p);

                } else if (arrayListSize ==4 && oxygenCounter>1) {
                    AtomType N_2p = new AtomType(Nitrogen.INSTANCE, "N_2p");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, N_2p);
                } else {
                    AtomType N_3p = new AtomType(Nitrogen.INSTANCE, "N_3p");
                    //System.out.println("The atom " + (i) +" is N_3 " );
                    atomIdentifierMap.put(i, N_3p);
                }

            } else if (retriveArrayFirstElementName.equals("S")) {
                //types are -S-, - -S=O =O, -S=O =O, -S=O =O -OH(-)
                int arrayListSize = connectivityModified.get(i).size();
                if(arrayListSize ==2){
                    AtomType S_2p = new AtomType(Sulfur.INSTANCE, "S_2p");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, S_2p);
                }else {
                    AtomType S_3p = new AtomType(Sulfur.INSTANCE, "S_3p");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, S_3p);
                }

            } else if (retriveArrayFirstElementName.equals("P")){
                AtomType P_3 = new AtomType(Phosphorus.INSTANCE, "P_3");
                //System.out.println("The atom " + (i+1) +" is P_3 " );
                atomIdentifierMap.put(i, P_3);

            } else if (retriveArrayFirstElementName.equals("SI")) {
                AtomType Si = new AtomType(Silicon.INSTANCE, "Si");
                // System.out.println("The atom " + (i+1) +" is Si " );
                atomIdentifierMap.put(i, Si);
                int arrayListSize = connectivityModified.get(i).size();

            }  else if (retriveArrayFirstElementName.equals("BR")) {
                AtomType Br = new AtomType(Bromine.INSTANCE, "Br");
                //System.out.println("The atom " + (i) +" is Br " );
                atomIdentifierMap.put(i, Br);

            } else if (retriveArrayFirstElementName.equals("F")) {
                AtomType F = new AtomType(Fluorine.INSTANCE, "F");
                // System.out.println("The atom " + (i+1) +" is F " );
                atomIdentifierMap.put(i, F);

            } else {//Metal ions
                if (retriveArrayFirstElementName.equals("RH")) {
                    AtomType Rh = new AtomType(Rhodium.INSTANCE, "Rh");
                    // System.out.println("The atom " + (i+1) +" is Rh " );
                    atomIdentifierMap.put(i, Rh);

                } else if (retriveArrayFirstElementName.equals("RU")) {
                    AtomType Ru = new AtomType(Ruthenium.INSTANCE, "Ru");
                    // System.out.println("The atom " + (i+1) +" is Ru " );
                    atomIdentifierMap.put(i, Ru);

                } else if (retriveArrayFirstElementName.equals("NI")) {
                    AtomType Ni = new AtomType(Nickel.INSTANCE, "Ni");
                    // System.out.println("The atom " + (i+1) +" is Ni " );
                    atomIdentifierMap.put(i, Ni);

                } else if (retriveArrayFirstElementName.equals("CU")) {
                    AtomType Cu = new AtomType(Copper.INSTANCE, "Cu");
                    //System.out.println("The atom " + (i+1) +" is Cu " );
                    atomIdentifierMap.put(i, Cu);

                } else if (retriveArrayFirstElementName.equals("FE")) { // Not this
                    int arrayListSize = connectivityModified.get(i).size();
                    if (arrayListSize == 3) {
                        // connected to two elements
                        AtomType Fe_2 = new AtomType(Iron.INSTANCE, "Fe_2");
                        // System.out.println("The atom " + (i+1) +" is Fe_2 " );
                        atomIdentifierMap.put(i, Fe_2);

                    } else {
                        AtomType Fe_3 = new AtomType(Iron.INSTANCE, "Fe_3");
                        // System.out.println("The atom " + (i+1) +" is N_2 " );
                        atomIdentifierMap.put(i, Fe_3);
                    }
                    //System.out.println(atomValues);

                } else if (retriveArrayFirstElementName.equals("CO")) {
                    AtomType Co = new AtomType(Cobalt.INSTANCE, "Co");
                    //System.out.println("The atom " + (i+1) +" is Co " );
                    atomIdentifierMap.put(i, Co);

                } else if (retriveArrayFirstElementName.equals("CR")) {
                    AtomType Cr = new AtomType(Chromium.INSTANCE, "Cr");
                    //System.out.println("The atom " + (i+1) +" is Cr " );
                    atomIdentifierMap.put(i, Cr);

                } else if (retriveArrayFirstElementName.equals("PD")) {
                    AtomType Pd = new AtomType(Palladium.INSTANCE, "Pd");
                    //System.out.println("The atom " + (i+1) +" is Pd " );
                    atomIdentifierMap.put(i, Pd);

                } else if (retriveArrayFirstElementName.equals("MO")) {
                    AtomType Mo = new AtomType(Molybdenum.INSTANCE, "Mo");
                    // System.out.println("The atom " + (i+1) +" is Mo " );
                    atomIdentifierMap.put(i, Mo);

                } else if (retriveArrayFirstElementName.equals("ZR")) {
                    AtomType Zr = new AtomType(Zirconium.INSTANCE, "Zr");
                    // System.out.println("The atom " + (i+1) +" is Zr " );
                    atomIdentifierMap.put(i, Zr);

                } else if (retriveArrayFirstElementName.equals("V")) {
                    AtomType V = new AtomType(Vanadium.INSTANCE, "V");
                    // System.out.println("The atom " + (i+1) +" is V " );
                    atomIdentifierMap.put(i, V);

                } else if (retriveArrayFirstElementName.equals("W")) {// Not this
                    AtomType W = new AtomType(Tungsten.INSTANCE, "W");
                    //System.out.println("The atom " + (i+1) +" is W " );
                    atomIdentifierMap.put(i, W);

                } else if (retriveArrayFirstElementName.equals("MG")) {
                    AtomType Mg = new AtomType(Magnesium.INSTANCE, "Mg");
                    // System.out.println("The atom " + (i+1) +" is Mg " );
                    atomIdentifierMap.put(i, Mg);

                } else if (retriveArrayFirstElementName.equals("ZN")) {
                    AtomType Zn = new AtomType(Zinc.INSTANCE, "Zn");
                    // System.out.println("The atom " + (i+1) +" is Zn " );
                    atomIdentifierMap.put(i, Zn);

                } else if (retriveArrayFirstElementName.equals("IR")) {
                    AtomType Ir = new AtomType(Iridium.INSTANCE, "Ir");
                    // System.out.println("The atom " + (i+1) +" is Zn " );
                    atomIdentifierMap.put(i, Ir);
                } else {
                    //System.out.println(retriveArrayFirstElementName);
                }
            }
        }
        //System.out.println(atomIdentifierMap + "atomIdentifier");
        return atomIdentifierMap;
    }
    public Map<String, AtomType> getTypeMapNew(){
        return typeMapNew;
    }
    public ArrayList<ArrayList<Integer>> getConnectivity(){
        //remove buildSpecies when running single file. When working along with
        // species= buildSpecies(confName);
        return connectivity;
    }

    public ArrayList<ArrayList<Integer>> getConnectivityWithSpecies(String confName){
        readPDBFile(confName);
        return connectivity;
    }
    public ArrayList<ArrayList<Integer>> getconnectivityModified (ArrayList<ArrayList<Integer>>connectivity){
        for (int i = 0; i < connectivity.size(); i++) {
            ArrayList<Integer> newList = new ArrayList<>();
            newList.add(i);
            for (int j = 1; j < connectivity.get(i).size(); j++) {
                newList.add(connectivity.get(i).get(j) - 1);
            }
            connectivityModified.add(newList);
        }
        return connectivityModified;
    }

    public Map<Integer,String> getAtomMap( ArrayList<ArrayList<Integer>> connectivity){
        return atomMap;
    }
    public HashMap<Integer, String> getatomMapModified( Map<Integer, String> atomMap){
        for (Map.Entry<Integer, String> entry : atomMap.entrySet()) {
            atomMapModified.put(entry.getKey() - 1, entry.getValue());
        }
        return atomMapModified;
    }
    public List<int[]> bondSorter(List<int[]> duplets,  Map<Integer, String> atomIdentifierMapModified){
        //List<int[]> priorityList = new ArrayList<>();
        List<int[]> dupletsActual = new ArrayList<>();
        for(int i=0; i<duplets.size(); i++){
            int[] newCombo = new int[2];
            int firstElement = duplets.get(i)[0];
            int secondElement = duplets.get(i)[1];
            String atomOne = String.valueOf(atomIdentifierMapModified.get(firstElement));
            String atomTwo = String.valueOf(atomIdentifierMapModified.get(secondElement));
            int numOne = priorityMapGenerator(atomOne);
            int numTwo = priorityMapGenerator(atomTwo);
            if(numOne == numTwo){
                newCombo = new int[]{firstElement, secondElement};
            } else if (numOne > numTwo) {
                newCombo = new int[]{secondElement, firstElement};
            } else {
                newCombo = new int[]{firstElement, secondElement};
            }
            dupletsActual.add(newCombo);
            int[] newBond = {numOne, numTwo};
            //priorityList.add(newBond);
        }
        //System.out.println(Arrays.deepToString(duplets.toArray()));
        //System.out.println(Arrays.deepToString(dupletsActual.toArray()));
        //System.out.println(atomIdentifierMapModified);
        //System.out.println(Arrays.deepToString(priorityList.toArray()));
        return dupletsActual;
    }
    public List<int[]> angleSorter(List<int[]> triplets,  Map<Integer, String> atomIdentifierMapModified){
        //List<int[]> priorityList = new ArrayList<>();
        List<int[]> tripletsActual = new ArrayList<>();
        for(int i=0; i<triplets.size(); i++){
            int[] newCombo = new int[2];
            String[] newComboString = new String[3];
            int firstElement = triplets.get(i)[0];
            int secondElement = triplets.get(i)[1];
            int thirdElement = triplets.get(i)[2];
            String atomOne = String.valueOf(atomIdentifierMapModified.get(firstElement));
            String atomTwo = String.valueOf(atomIdentifierMapModified.get(secondElement));
            String atomThree = String.valueOf(atomIdentifierMapModified.get(thirdElement));
            int numOne = priorityMapGenerator(atomOne);
            int numThree = priorityMapGenerator(atomThree);
            if (numOne > numThree) {
                // newComboString = new String[]{atomThree, atomTwo, atomOne};
                newCombo = new int[]{thirdElement,secondElement, firstElement};
                // System.out.println(Arrays.toString(newCombo) + " part1 " + Arrays.toString(newComboString));
            } else {
                // newComboString = new String[]{atomThree, atomTwo, atomOne};
                newCombo = new int[]{firstElement,secondElement, thirdElement};
                // System.out.println(Arrays.toString(newCombo) + " part2 " + Arrays.toString(newComboString));
            }
            tripletsActual.add(newCombo);
            //int[] newBond = {numOne, numThree};
            //priorityList.add(newBond);
        }
        System.out.println(Arrays.deepToString(triplets.toArray()));
        System.out.println(Arrays.deepToString(tripletsActual.toArray()));
        System.out.println(atomIdentifierMapModified);
        //System.out.println(Arrays.deepToString(priorityList.toArray()));
        return tripletsActual;
    }
    public List<int[]> torsionSorter(List<int[]> quadruplets,  Map<Integer, String> atomIdentifierMapModified){
        List<int[]> quadrupletsActual = new ArrayList<>();
        for(int i=0; i<quadruplets.size(); i++){
            int[] newCombo = new int[4];
            int firstElement = quadruplets.get(i)[0];
            int secondElement = quadruplets.get(i)[1];
            int thirdElement = quadruplets.get(i)[2];
            int fourthElement = quadruplets.get(i)[3];
            String atomOne = String.valueOf(atomIdentifierMapModified.get(firstElement));
            String atomFour = String.valueOf(atomIdentifierMapModified.get(fourthElement));
            int numOne = priorityMapGenerator(atomOne);
            int numFour = priorityMapGenerator(atomFour);
            if(numOne == numFour){
                newCombo = new int[]{firstElement, secondElement, thirdElement, fourthElement};
            } else if (numOne > numFour) {
                newCombo = new int[]{fourthElement, thirdElement,secondElement, firstElement};
            } else {
                newCombo = new int[]{firstElement, secondElement, thirdElement, fourthElement};
            }
            quadrupletsActual.add(newCombo);
        }
        System.out.println(Arrays.deepToString(quadruplets.toArray()));
        System.out.println(Arrays.deepToString(quadrupletsActual.toArray()));
        System.out.println(atomIdentifierMapModified);
        return quadrupletsActual;
    }

    public Map<String[],List<int[]>> idenBondTypes(List<int[]> duplets, Map<Integer, String> atomIdentifierMapModified){
        System.out.println(Arrays.deepToString(duplets.toArray()) + " idenBond");
        Map<String[],List<int[]>> bondTypesMap = new HashMap<>();
        // System.out.println(atomIdentifierMapModified);
        ArrayList<String[]> bondTypes = new ArrayList<>();
        ArrayList<Integer> bondValues = new ArrayList<>();
        //System.out.println(bondValues);
        //System.out.println("\nThe printing of the bond types");
        //System.out.println("The printing of the bond types");
        //form bonded pairs of atoms
        for(int i =0; i< duplets.size(); i++){
            int firstElement = duplets.get(i)[0];
            int secondElement = duplets.get(i)[1];
            int bondValue = bondsNum.get(i);
            // System.out.println(bondValue);
            String atomOne = String.valueOf(atomIdentifierMapModified.get(firstElement));
            String atomTwo = String.valueOf(atomIdentifierMapModified.get(secondElement));
            String[] newBond = {atomOne, atomTwo};
            //Arrays.sort(newBond); // Sort the elements
            boolean isPresent = false;
            for (String[] array : bondTypes) {
                if (array.length == newBond.length) {
                    boolean isEqual = true;
                    for (int j = 0; j < array.length; j++) {
                        if (!array[j].equals(newBond[j])) {
                            isEqual = false;
                            //System.out.println("not equal " + newBond[j]);
                            break;

                        }
                    }
                    if (isEqual) {
                        isPresent = true;
                        // System.out.println("is equal");
                        break;
                    }
                }
            }
            // If target array is not present, add it to the ArrayList
            if (!isPresent) {
                bondValues.add(bondValue);
                bondTypes.add(newBond);
            }

        }


        // System.out.println(Arrays.deepToString(bondTypes.toArray()) + " bondTypes");
        // System.out.println(bondValues + " bonds");
        //identify the pairs of individual types
        for(int i=0; i<bondTypes.size();i++){
            System.out.println(bondTypes.get(i)[0] + " " + bondTypes.get(i)[1]);
            List<int[]> dupletsNew = new ArrayList<>();
            for(int k =0; k< duplets.size();k++){
                int firstElement = duplets.get(k)[0];
                int secondElement = duplets.get(k)[1];
                String atomOne = String.valueOf(atomIdentifierMapModified.get(firstElement));
                String atomTwo = String.valueOf(atomIdentifierMapModified.get(secondElement));
                //System.out.println(firstElement + " " + atomOne +" "+ secondElement + " " + atomTwo);
                if((bondTypes.get(i)[0].equals(atomOne) && bondTypes.get(i)[1].equals(atomTwo))||(bondTypes.get(i)[0].equals(atomTwo) && bondTypes.get(i)[1].equals(atomOne)) ){
                    int[] sameBond = {firstElement, secondElement};
                    dupletsNew.add(sameBond);
                }
                // System.out.println(Arrays.deepToString(dupletsNew.toArray()) + " dupletsNew After");
            }
            if(bondTypesMap.containsKey(bondTypes.get(i))){
                break;
            } else {
                bondTypesMap.put(bondTypes.get(i), dupletsNew);
            }
        }
        for (Map.Entry<String[], List<int[]>> entry : bondTypesMap.entrySet()) {
            String[] bondType = entry.getKey();
            List<int[]> bonds = entry.getValue();
            System.out.println(Arrays.toString(bondType) + ": " + Arrays.deepToString(bonds.toArray()));
        }
        System.out.println(bondValues);
        System.out.println("\n");
        return bondTypesMap;
    }

    public Map<String, Integer> getangleMap(List<int[]>dupletsSorted, ArrayList<Integer> bondsNum){
        Map<String, Integer> anglemap = new HashMap<>();
        Map<String, Integer> map = new HashMap<>();

        // Iterate through the List<int[]> and ArrayList<Integer>
        // to create the mapping
        for (int i = 0; i < dupletsSorted.size(); i++) {
            String key = Arrays.toString(dupletsSorted.get(i));
            Integer value = bondsNum.get(i);
            map.put(key, value);
        }
        int[] arr = {0, 1};
        System.out.println(anglemap.get(arr));
        return anglemap;
    }


    public Map<String[],List<int[]>> idenAngleTypes( List<int[]> triplets,  Map<Integer, String> atomIdentifierMapModified){
        // System.out.println(Arrays.deepToString(triplets.toArray()));
        Map<String[],List<int[]>> angleTypesMap = new HashMap<>();
        //System.out.println(atomIdentifierMapModified);
        ArrayList<String[]> angleTypes = new ArrayList<>();
        System.out.println("The printing of the angle types");
        //form bonded pairs of atoms
        for(int i =0; i< triplets.size(); i++){
            int firstElement = triplets.get(i)[0];
            int secondElement = triplets.get(i)[1];
            int thirdElement = triplets.get(i)[2];
            String atomOne = String.valueOf(atomIdentifierMapModified.get(firstElement));
            String atomTwo = String.valueOf(atomIdentifierMapModified.get(secondElement));
            String atomThree = String.valueOf(atomIdentifierMapModified.get(thirdElement));
            String[] newAngle = {atomOne, atomTwo, atomThree};
            //Arrays.sort(newAngle);
            boolean isPresent = false;
            for (String[] array : angleTypes) {
                //System.out.println(Arrays.toString(array));
                if (array.length == newAngle.length) {
                    boolean isEqual = true;
                    for (int j = 0; j < array.length; j++) {
                        if (!array[j].equals(newAngle[j])) {
                            isEqual = false;
                            //System.out.println("not equal " + newAngle[j]);
                            break;
                        }
                    }
                    if (isEqual) {
                        isPresent = true;
                        //System.out.println("is equal");
                        break;
                    }
                }
            }
            // If target array is not present, add it to the ArrayList
            if (!isPresent) {
                angleTypes.add(newAngle);
            }
        }
        System.out.println(Arrays.deepToString(angleTypes.toArray()) + " angleTypes");
        //identify the pairs of individual types
        for(int i=0; i<angleTypes.size();i++){
            //System.out.println(angleTypes.get(i)[0] + " " + angleTypes.get(i)[1]+ " " + angleTypes.get(i)[2]);
            List<int[]> tripletNew = new ArrayList<>();
            for(int k =0; k< triplets.size();k++){
                int firstElement = triplets.get(k)[0];
                int secondElement = triplets.get(k)[1];
                int thirdElement = triplets.get(k)[2];
                String atomOne = String.valueOf(atomIdentifierMapModified.get(firstElement));
                String atomTwo = String.valueOf(atomIdentifierMapModified.get(secondElement));
                String atomThree = String.valueOf(atomIdentifierMapModified.get(thirdElement));
                // System.out.println(firstElement + " " + atomOneActual +" "+ secondElement + " " + atomTwoActual);
                if((angleTypes.get(i)[0].equals(atomOne) && angleTypes.get(i)[1].equals(atomTwo)&& angleTypes.get(i)[2].equals(atomThree))||(angleTypes.get(i)[0].equals(atomThree) && angleTypes.get(i)[1].equals(atomTwo)&& angleTypes.get(i)[2].equals(atomOne))){
                    int[] sameAngle = {firstElement, secondElement, thirdElement};
                    tripletNew.add(sameAngle);
                }
                //System.out.println(Arrays.deepToString(dupletsNew.toArray()) + " dupletsNew After");
            }
            if(angleTypesMap.containsKey(angleTypes.get(i))){
                break;
            } else {
                angleTypesMap.put(angleTypes.get(i), tripletNew);
            }
        }
        for (Map.Entry<String[], List<int[]>> entry : angleTypesMap.entrySet()) {
            String[] bondType = entry.getKey();
            List<int[]> angles = entry.getValue();
            System.out.println(Arrays.toString(bondType) + ": " + Arrays.deepToString(angles.toArray()));
        }
        System.out.println("\n");
        return angleTypesMap;
    }

    public Map<String[],List<int[]>> idenTorsionTypes( List<int[]> quadruplets,  Map<Integer, String> atomIdentifierMapModified){
        //System.out.println(Arrays.deepToString(quadruplets.toArray()) + "quadruplets");
        Map<String[],List<int[]>> torsionTypesMap = new HashMap<>();
        //System.out.println(atomIdentifierMapModified);
        ArrayList<String[]> torsionTypes = new ArrayList<>();
        //System.out.println("The printing of the torsion types");
        //form bonded pairs of atoms
        for(int i =0; i< quadruplets.size(); i++){
            int firstElement = quadruplets.get(i)[0];
            int secondElement = quadruplets.get(i)[1];
            int thirdElement = quadruplets.get(i)[2];
            int fourthElement = quadruplets.get(i)[3];
            String atomOne = String.valueOf(atomIdentifierMapModified.get(firstElement));
            String atomTwo = String.valueOf(atomIdentifierMapModified.get(secondElement));
            String atomThree = String.valueOf(atomIdentifierMapModified.get(thirdElement));
            String atomFour = String.valueOf(atomIdentifierMapModified.get(fourthElement));
            String[] newTorsion = {atomOne, atomTwo, atomThree, atomFour};
            boolean isPresent = false;
            for (String[] array : torsionTypes) {
                if (array.length == newTorsion.length) {
                    boolean isEqual = true;
                    for (int j = 0; j < array.length; j++) {
                        if (!array[j].equals(newTorsion[j])) {
                            isEqual = false;
                            //System.out.println("not equal " + newAngle[j]);
                            break;
                        }
                    }
                    if (isEqual) {
                        isPresent = true;
                        //System.out.println("is equal");
                        break;
                    }
                }
            }
            // If target array is not present, add it to the ArrayList
            if (!isPresent) {
                torsionTypes.add(newTorsion);
            }
        }
        // System.out.println(Arrays.deepToString(torsionTypes.toArray()) + " torsionTypes");
        //identify the pairs of individual types
        for(int i=0; i<torsionTypes.size();i++){
            //System.out.println(torsionTypes.get(i)[0] + " " + torsionTypes.get(i)[1]+ " " + torsionTypes.get(i)[2]);
            List<int[]> quadrupletNew = new ArrayList<>();
            for(int k =0; k< quadruplets.size();k++){
                // System.out.println(Arrays.toString(quadruplets.get(k)) + " quadruplets");
                int firstElement = quadruplets.get(k)[0];
                int secondElement = quadruplets.get(k)[1];
                int thirdElement = quadruplets.get(k)[2];
                int fourthElement = quadruplets.get(k)[3];
                String atomOne = String.valueOf(atomIdentifierMapModified.get(firstElement));
                String atomTwo = String.valueOf(atomIdentifierMapModified.get(secondElement));
                String atomThree = String.valueOf(atomIdentifierMapModified.get(thirdElement));
                String atomFour = String.valueOf(atomIdentifierMapModified.get(fourthElement));
                // System.out.println(firstElement + " " + atomOneActual +" "+ secondElement + " " + atomTwoActual);
                if((torsionTypes.get(i)[0].equals(atomOne) && torsionTypes.get(i)[1].equals(atomTwo)&& torsionTypes.get(i)[2].equals(atomThree)&& torsionTypes.get(i)[3].equals(atomFour))||(torsionTypes.get(i)[0].equals(atomFour) && torsionTypes.get(i)[1].equals(atomThree)&& torsionTypes.get(i)[2].equals(atomTwo))&& torsionTypes.get(i)[3].equals(atomOne)){
                    int[] sameAngle = {firstElement, secondElement, thirdElement, fourthElement};
                    quadrupletNew.add(sameAngle);
                }
                //System.out.println(Arrays.deepToString(dupletsNew.toArray()) + " dupletsNew After");
            }
            if(torsionTypesMap.containsKey(torsionTypes.get(i))){
                break;
            } else {
                torsionTypesMap.put(torsionTypes.get(i), quadrupletNew);
            }
        }
        for (Map.Entry<String[], List<int[]>> entry : torsionTypesMap.entrySet()) {
            String[] torsionType = entry.getKey();
            List<int[]> torsions = entry.getValue();
            System.out.println(Arrays.toString(torsionType) + ": " + Arrays.deepToString(torsions.toArray()));
        }
        System.out.println("\n");
        return torsionTypesMap;
    }

    public List<int[]> getBondedAtomList (ArrayList<ArrayList<Integer>> connectivityModified){
        ArrayList<ArrayList<Integer>> listOfBonds = new ArrayList<>();
        for (ArrayList<Integer> list : connectivityModified) {
            int firstElement = list.get(0);
            for (int i = 1; i < list.size(); i++) {
                ArrayList<Integer> arrayList = new ArrayList<>();
                int secondElement = list.get(i);
                if(firstElement<secondElement){
                    arrayList.add(firstElement);
                    arrayList.add(list.get(i));
                    listOfBonds.add(arrayList);
                    //System.out.println(listOfBonds + " listOfBonds");
                }
            }
        }
        for (List<Integer> innerList : listOfBonds) {
            int[] array = innerList.stream().mapToInt(Integer::intValue).toArray();
            intListArray.add(array);
        }
        return intListArray;
    }
    public ArrayList<ArrayList<ArrayList<Integer>>> generateCombinationsForNestedList(ArrayList<ArrayList<Integer>> inputList) {
        ArrayList<ArrayList<ArrayList<Integer>>> outputList = new ArrayList<>();
        for (ArrayList<Integer> innerList : inputList) {
            ArrayList<ArrayList<Integer>> innerOutputList = generateCombinations(innerList);
            outputList.add(innerOutputList);
        }
        return outputList;
    }
    public ArrayList<ArrayList<Integer>> generateCombinations(List<Integer> inputList) {
        ArrayList<ArrayList<Integer>> outputList = new ArrayList<>();
        if (inputList.size() < 3) {
            return outputList;
        }
        int center = inputList.get(0);
        for (int i = 1; i < inputList.size(); i++) {
            for (int j = i + 1; j < inputList.size(); j++) {
                ArrayList<Integer> tempList = new ArrayList<>(Arrays.asList(inputList.get(i), center, inputList.get(j)));
                outputList.add(tempList);
            }
        }
        return outputList;
    }

    public List<int[]> getAngleList(ArrayList<ArrayList<Integer>> connectivityModified) {
        ArrayList<ArrayList<ArrayList<Integer>>> arrayListOfAngles = generateCombinationsForNestedList(connectivityModified);


        for (ArrayList<ArrayList<Integer>> subList1 : arrayListOfAngles) {
            for (ArrayList<Integer> subList2 : subList1) {
                int[] arr = subList2.stream().mapToInt(Integer::intValue).toArray();
                tripletsList.add(arr);
            }
        }
        return tripletsList ;
    }
    public ArrayList<Integer> getlistOfTorsions(ArrayList<ArrayList<Integer>> connectivityModified){
        ArrayList<Integer> torsionPairs = new ArrayList<>();
        ArrayList<Integer> temp = new ArrayList<>();
        for (int i =0; i<connectivityModified.size(); i++){
            ArrayList<Integer> connectivityElement;
            connectivityElement = connectivityModified.get(i);
            // System.out.println(connectivityElement);
            if(connectivityElement.size() > 2){
                for (int j =1; j <connectivityElement.size(); j++){
                    int individualValue = connectivityElement.get(j);
                    ArrayList<Integer> connectedElement = new ArrayList<>();
                    connectedElement = connectivityModified.get(individualValue-1);
                    if(connectivityElement.get(0) < connectivityElement.get(j)){
                        if(connectedElement.size() > 2 ){
                            torsionPairs.add(connectivityElement.get(0));
                            torsionPairs.add(connectedElement.get(0));
                            // System.out.println(torsionPairs);
                        }
                    }
                }
            }
        }
        //System.out.println(torsionPairs + "Torsionpairs");
        return torsionPairs;
    }
    public List<List<Integer>> setTorsionPairs(ArrayList<ArrayList<Integer>> connectivityModified){
        List<List<Integer>> listOfTorsions = new ArrayList<>();
        ArrayList<Integer> torsionPairs = getlistOfTorsions(connectivityModified);
        for(int i = 0; i < torsionPairs.size(); i++){
            int elementOne = torsionPairs.get(i);
            int elementTwo = torsionPairs.get(i+1);
            List<Integer> mainElement = connectivityModified.get(elementOne-1);
            List<Integer> secondElement = connectivityModified.get(elementTwo-1);
            i++;
            for(int j=1; j<mainElement.size(); j++){
                int initialElement = mainElement.get(j);
                if(initialElement != elementOne && initialElement != elementTwo){
                    for(int k =1; k<secondElement.size(); k++){
                        List<Integer> temp = new ArrayList<>();
                        int finalElement = secondElement.get(k);
                        if(finalElement != elementOne && finalElement != elementTwo){
                            temp.add(initialElement);
                            temp.add(elementOne);
                            temp.add(elementTwo);
                            temp.add(finalElement);
                            listOfTorsions.add(temp);
                        }
                    }
                }
            }
        }
        //System.out.println(listOfTorsions);
        return listOfTorsions;
    }

    public List<int[]> getTorsionList(ArrayList<ArrayList<Integer>> connectivity) {
        List<List<Integer>> listOfTorsion=setTorsionPairs(connectivity);
        quadruplets = listOfTorsion.stream()
                .map(sublist -> new int[]{sublist.get(0) - 1, sublist.get(1) - 1, sublist.get(2) - 1, sublist.get(3) - 1})
                .collect(Collectors.toList());
        //System.out.println(quadruplets);
        return quadruplets;

    }


    public Set<String> uniqueElementIdentifier(){
        Set<String> uniqueAtoms = new HashSet<>();
        for(int i =0; i<atomIdentifierMap.size(); i++){
            AtomType atomName = atomIdentifierMap.get(i);
            uniqueAtoms.add(String.valueOf(atomName));
        }
        return uniqueAtoms;
    }

    public Map<String, double[]> atomicPotMap(){
        Map<String, double[]> atomicPot = new HashMap<>();
        Set<String> uniqueElement= uniqueElementIdentifier();
        Iterator<String> iterator = uniqueElement.iterator();
        while (iterator.hasNext()){
            String str = iterator.next();
            //System.out.println(str + "str");
            String elementName = str.substring(str.indexOf("[") + 1, str.indexOf("]"));
            System.out.println( str+ " " +elementName);
            double[] atomicPotValues = atomicPot(elementName);
            //System.out.println(str + " " + Arrays.toString(atomicPotValues));
            //atomicPotMap.put(str, atomicPotValues);
            atomicPot.put(elementName, atomicPotValues);
            System.out.println(elementName + " " + Arrays.toString(atomicPotValues));
        }
        return atomicPot;
    }

    public Map<Integer, Integer> coordinationNumberDeterminer(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomIdentifierMapModified){

        Map<Integer, Integer> coordinationMap = new HashMap<>();
        System.out.println(connectivityModified);
        System.out.println(atomIdentifierMapModified);
        for(int i =0; i<connectivityModified.size();i++){
            int size = connectivityModified.get(i).size() -1;
            //System.out.println(i+" "+ connectivityModified.get(i)+" "+ size+" "+atomIdentifierMapModified.get(i));
            if(atomIdentifierMapModified.get(i).equals("C_2") || atomIdentifierMapModified.get(i).equals("C_1")){
                coordinationMapOutput.put(i,3);
            } else {
                coordinationMapOutput.put(i,0);
            }
            coordinationMap.put(i, size);
        }
        // System.out.println(coordinationMap+ " coordinationMap");
        // System.out.println(coordinationMapOutput+ " coordinationMap");
        return coordinationMapOutput;
    }

    public double [] atomicPot (String atomtype){
        HashMap<String, double[]> atomicConstant = new HashMap<>();
        atomicConstant.put("V", new double[]{1.402, 109.47, 3.144, 0.016, 12.0, 2.679,3.65,0.7});
        atomicConstant.put("Co", new double[]{1.241, 90.0, 2.872, 0.014, 12.0, 2.43,4.105,0.7});
        atomicConstant.put("Mo", new double[]{1.467,90,3.052,0.056,12.0,3.40,3.465,0.2});
        atomicConstant.put("Fe_2", new double[]{1.27, 109.47, 2.912, 0.013, 12.0, 2.43,0.7});
        atomicConstant.put("Fe_3", new double[]{1.335, 90.0, 2.912, 0.013, 12.0, 2.43,3.76,0.7});
        atomicConstant.put("Mn", new double[]{1.382, 90.0, 2.961, 0.013, 12.0, 2.43,3.325,0.7});
        atomicConstant.put("Ni", new double[]{1.164, 90.0, 2.834, 0.015, 12.0, 2.43,4.465,0.7});
        atomicConstant.put("Pd", new double[]{1.338, 90.0, 2.899, 0.048, 12.0, 3.21,4.32,0.2});
        atomicConstant.put("Zr", new double[]{1.564, 109.47, 3.124, 0.069, 12.0, 3.667,3.40,0.2});
        atomicConstant.put("W+6", new double[]{1.392, 90.0, 3.069, 0.067, 12.0, 3.7,4.63,0.1});
        atomicConstant.put("W+4", new double[]{1.526, 109.47, 3.069, 0.067, 12.0, 3.7,4.63,0.1});
        atomicConstant.put("Zn", new double[]{1.193, 109.47, 2.763, 0.124, 12.0, 1.308,5.106,0.7});
        atomicConstant.put("Cu", new double[]{1.302, 109.47, 3.495, 0.005, 12.0, 1.756,3.729,0.7});
        atomicConstant.put("Rh", new double[]{1.332, 90.0, 2.929, 0.053, 12.0, 3.5,3.975,0.2});
        atomicConstant.put("Cr", new double[]{1.345, 90.0, 3.023, 0.015, 12.0, 2.463,3.415,0.7});
        atomicConstant.put("Ru", new double[]{1.478, 90.0, 2.963, 0.056, 12.0, 3.4,3.575,0.2});
        atomicConstant.put("Si", new double[]{1.117, 109.47, 4.295, 0.402, 12.175, 2.323,4.168,1.225});
        atomicConstant.put("C_Ar", new double[]{0.729, 120.0, 3.851, 0.105, 12.73, 1.912, 5.343, 2.0}); //C_Ar
        atomicConstant.put("C_3", new double[]{0.757, 109.47, 3.851, 0.105, 12.73, 1.912, 5.343,2.0});
        atomicConstant.put("C_2", new double[]{0.732, 120.0, 3.851, 0.105, 12.73, 1.912, 5.343,2.0});
        atomicConstant.put("C_1", new double[]{0.706, 180, 3.851,0.105,12.73, 1.912, 5.343,2.0});
        atomicConstant.put("N_2", new double[]{0.699, 120.0, 3.66, 0.069, 13.407, 2.544,6.899,2.0});
        atomicConstant.put("N_3", new double[]{0.7, 106.7, 3.66, 0.069, 13.407, 2.544, 6.899,0.45});
        atomicConstant.put("N_1", new double[]{0.685, 111.2, 3.66, 0.069, 13.407, 2.544,6.899,2.0});
        atomicConstant.put("O_Ar", new double[]{0.68, 110.0, 3.5, 0.06, 14.085, 2.3, 8.741, 2.0});
        atomicConstant.put("O_2", new double[]{0.68, 110.0, 3.5, 0.06, 14.085, 2.3, 8.741, 2.0});
        atomicConstant.put("O_3", new double[]{0.658, 104.51, 3.5, 0.06, 14.085, 2.3,8.74,0.018});
        atomicConstant.put("O_1", new double[]{0.639, 180.0, 3.5, 0.06, 14.085, 2.3, 8.741,2.0});
        atomicConstant.put("S_2", new double[]{1.077, 92.2, 4.035, 0.274, 13.969, 2.703,6.928,1.25});
        atomicConstant.put("S_3", new double[]{1.064, 92.1, 4.035, 0.274, 13.969, 2.703,6.928,0.484});
        atomicConstant.put("C_Arp", new double[]{0.729, 120.0, 3.851, 0.105, 12.73, 1.912, 5.343, 2.0}); //C_Ar
        atomicConstant.put("C_3p", new double[]{0.757, 109.47, 3.851, 0.105, 12.73, 1.912, 5.343,2.0});
        atomicConstant.put("C_2p", new double[]{0.732, 120.0, 3.851, 0.105, 12.73, 1.912, 5.343,2.0});
        atomicConstant.put("C_1p", new double[]{0.706, 180, 3.851,0.105,12.73, 1.912, 5.343,2.0});
        atomicConstant.put("N_2p", new double[]{0.699, 120.0, 3.66, 0.069, 13.407, 2.544,6.899,2.0});
        atomicConstant.put("N_3p", new double[]{0.7, 106.7, 3.66, 0.069, 13.407, 2.544, 6.899,0.45});
        atomicConstant.put("N_1p", new double[]{0.685, 111.2, 3.66, 0.069, 13.407, 2.544,6.899,2.0});
        atomicConstant.put("O_Arp", new double[]{0.68, 110.0, 3.5, 0.06, 14.085, 2.3, 8.741, 2.0});
        atomicConstant.put("O_2p", new double[]{0.68, 110.0, 3.5, 0.06, 14.085, 2.3, 8.741, 2.0});
        atomicConstant.put("O_3p", new double[]{0.658, 104.51, 3.5, 0.06, 14.085, 2.3,8.74,0.018});
        atomicConstant.put("O_1p", new double[]{0.639, 180.0, 3.5, 0.06, 14.085, 2.3, 8.741,2.0});
        atomicConstant.put("S_2p", new double[]{1.077, 92.2, 4.035, 0.274, 13.969, 2.703,6.928,1.25});
        atomicConstant.put("S_3p", new double[]{1.064, 92.1, 4.035, 0.274, 13.969, 2.703,6.928,0.484});
        atomicConstant.put("S+4", new double[]{1.049, 103.2, 4.035, 0.274, 13.969, 2.703,6.928,0.484});
        atomicConstant.put("S+6", new double[]{1.027, 109.47, 4.035, 0.274, 13.969, 2.703,6.928,0.484});
        atomicConstant.put("P_3", new double[]{1.101, 93.8, 4.147, 0.305, 13.072, 2.863,5.463,2.4});
        //atomicConstant.put("P+5", new double[]{1.056, 109.47, 4.147, 0.305, 13.072, 2.863,5.463});
        atomicConstant.put("H", new double[]{0.354, 180.0, 2.886, 0.044, 12.0, 0.71,4.528, 0.0});
        atomicConstant.put("H_p", new double[]{0.354, 180.0, 2.886, 0.044, 12.0, 0.71,4.528, 0.0});
        atomicConstant.put("Cl", new double[]{1.044,180.0,3.947,0.227,14.866,2.348,8.564,1.25});
        atomicConstant.put("K", new double[]{1.953, 180, 3.812, 0.035, 12, 1.165,2.421,0.7});
        atomicConstant.put("Br", new double[]{1.192, 180, 4.189, 0.251, 15, 2.519,7.790,0.7});
        atomicConstant.put("Na", new double[]{1.539, 180, 2.983, 0.03, 12.0 , 1.081,2.843,1.25});
        atomicConstant.put("Mg", new double[]{1.421,109.47, 3.021, 0.111, 12.0, 1.787,1.25});
        atomicConstant.put("Al", new double[]{1.244, 109.47, 4.499, 0.505, 11.278, 1.792,4.06,1.25});
        atomicConstant.put("Ca", new double[]{1.761, 90,3.399,0.238,12,2.141,3.231,0.7});
        atomicConstant.put("Ti3", new double[]{1.412,109.47,3.175,0.017,12.0,2.659,3.47,0.7});
        atomicConstant.put("Ti6", new double[]{1.412,90,3.175,0.017,12.0,2.659,0});
        atomicConstant.put("I", new double[]{1.382,180,4.5,0.339,15,2.65,6.822,0.2});
        atomicConstant.put("F", new double[]{0.668,180,3.364,0.050,14.762,1.735,10.874,2.0});
        atomicConstant.put("He", new double[]{0.849, 90.0, 2.362,0.056,15.24,0.098});
        atomicConstant.put("Ne", new double[]{0.920, 90.0, 3.243, 0.042, 15.440, 0.194});
        atomicConstant.put("Ar", new double[]{1.032, 90.0, 3.868,0.185,15.763, 0.300});
        double [] sample = atomicConstant.get(atomtype);
        return sample;
    }

    public double [] electronicParam (String atomtype){
        HashMap<String, double[]> electronicConstant = new HashMap<>();
        electronicConstant.put("C", new double[]{5.343,10.126,0.759,0.8563});
        electronicConstant.put("O", new double[]{8.741,13.364,0.669,0.9745});
        electronicConstant.put("N", new double[]{6.8991,1.760,0.715,0.9089});
        electronicConstant.put("F", new double[]{10.874,14.948,0.706,0.9206});
        electronicConstant.put("Na", new double[]{2.843,4.592,2.085,0.4364});
        electronicConstant.put("Si", new double[]{4.168,6.974,1.176,0.7737});
        electronicConstant.put("P", new double[]{5.463,8.000,1.102,0.8257});
        electronicConstant.put("S", new double[]{6.928,8.972,1.047,0.8690});
        electronicConstant.put("Cl", new double[]{8.564,9.892,0.994,0.9154});
        electronicConstant.put("K", new double[]{2.421,3.84,2.586,0.4524});
        electronicConstant.put("Br", new double[]{7.790,8.850,1.141,1.0253});
        electronicConstant.put("I", new double[]{6.822,7.524,1.333,1.0726});
        electronicConstant.put("H", new double[]{4.5280,13.8904,0.371,1.0698});
        double [] sample = electronicConstant.get(atomtype);
        return sample;
    }

    public ArrayList<ArrayList<Integer>> getConnectivityWithoutRunning(){
        return connectivity;
    }

    public ArrayList<ArrayList<Integer>> getConnectivityModifiedWithoutRunning(){
        return connectivityModified;
    }
    public Map<Integer,String> getAtomMapWithoutRunning(){
        return atomMap;
    }

    public HashMap<Integer, String> getAtomMapModifiedWithoutRunning(){
        return atomMapModified;
    }

    public List<int[]> getduplets(){return intListArray;}
    public List<int[]> gettriplets(){return tripletsList;}
    public List<int[]> getquadruplets(){return quadruplets;}
    public Map<Integer, String> getModifiedAtomIdentifierMap(){
        System.out.println(modifiedAtomIdentifierMap);
        return  modifiedAtomIdentifierMap;
    }
    public void main(String[] args) {
        String confName = "F://Avagadro//mopstrands//metallic//co4c24cl4";
    }
}
