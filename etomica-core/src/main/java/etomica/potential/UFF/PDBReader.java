package etomica.potential.UFF;


import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.molecule.IMolecule;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.util.collections.IntArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class PDBReader {
    static double molecularWeight =0;
    protected static Map<Integer,Vector> positions = new HashMap<>();
    protected static ISpecies sm;
    protected static Map<Integer, String> atomMap = new HashMap<>();
    public static ArrayList<ArrayList<Integer>> connectivity = new ArrayList<>();
    //static List<List> locationList = new ArrayList<>();
    public static HashMap<Integer, String> atomMapModified = new HashMap<>();
    public static Map<String, AtomType> typeMapNew = new HashMap<>();
   public static HashMap<Integer, String> atomMapNew = new HashMap<>();
    public static ArrayList<ArrayList<Integer>> connectivityModified = new ArrayList<>();
    public static ArrayList<ArrayList<Integer>> listOfBonds = new ArrayList<>();
    public static List<List<Integer>> listOfTorsions = new ArrayList<>();
    public static ArrayList<Integer> bondList = new ArrayList<>();
    public static Map<Integer, AtomType> atomIdentifierMap = new HashMap<>();
    static ISpecies species, speciesNew;

    public static void buildSpecies(String confName) {
        SpeciesBuilder speciesBuilder =  new SpeciesBuilder(Space3D.getInstance());
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
            String line = "";

            while ((line = bufReader.readLine()) != null) {

                parseLineReader(line, typeMap, atomMap, positions);

                if (line.startsWith("CONECT")) {
                    String[] parts = line.trim().split("\\s+");
                    int atomNumber = Integer.parseInt(parts[1]);
                    if (currentAtomList == null || atomNumber != currentAtomList.get(0)) {
                        // start a new list for the current atom
                        currentAtomList = new ArrayList<Integer>();
                        connectivity.add(currentAtomList);
                        currentAtomList.add(atomNumber);
                        //System.out.println(currentAtomList + "up");
                    }
                    for (int i = 2; i < parts.length; i++) {
                        currentAtomList.add(Integer.parseInt(parts[i]));
                        //System.out.println(currentAtomList + "down");
                    }
                }
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        String  molecularWeightString = String.format("%.4f", molecularWeight);
        System.out.println( "Molecular Weight: "+molecularWeightString);
        System.out.println(typeMap + "typeMap");
        System.out.println(connectivity + " connectivity");
        System.out.println(atomMap+" atomMap");
        System.out.println("Building done");
    }

    protected static void parseLineReader(String line, Map<String, AtomType > typeMap, Map<Integer, String> atomMap , Map<Integer, Vector> positions) {
        line = line.trim();
        int currentValue = -1;
        if (line.length() < 6) {
            return;
        }
        if (line.substring(0, 6).equals("HETATM")) {
            //coordinates of atom and create atom
            double x = Double.parseDouble(line.substring(30, 38));
            double y = Double.parseDouble(line.substring(38, 46));
            double z = Double.parseDouble(line.substring(46, 54));
            Vector positn = Vector.of(x, y, z);
            String symbol = line.substring(12, 16).trim();
            int atomNumber = Integer.parseInt(line.substring(8,11).trim());
            System.out.println(atomNumber + " " + positn);
            positions.put(atomNumber, positn);
            AtomType type;
            String symbolNew = line.substring(76, 78).trim();
            double weight = getMolecularWeight(symbolNew);

            // Add the atomic weight to the molecular weigh
            molecularWeight += weight;
            if (typeMap.containsKey(symbol)) {
                type = typeMap.get(symbol);
                //System.out.println(type + "contains");
            } else {
                type = new AtomType(new ElementSimple(line.substring(76, 78).trim()));
                //System.out.println(type + " newtype");
                typeMap.put(symbol, type);
            }
            if (atomMap.containsKey(atomNumber)){
                // System.out.println("Error");
            } else {
                atomMap.put(atomNumber, symbol);
            }
        }
    }

    public static <IAtomType> ISpecies getSpeciesNew (String confName){
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        AtomType typeNew;
        buildSpecies(confName);
       ArrayList<Integer> bondList = PDBReader.setBondList(connectivity, atomMap);
        System.out.println(bondList + " bondList");
        System.out.println(atomMap+ "  " + atomMap.size());
        Set<String> uniqueElements = PDBReader.uniqueElementIdentifierWithInputs(connectivity, atomMap);
        System.out.println(uniqueElements + "uniqueElements");
        Map<Integer, String> atomIdentifierMapModified = PDBReader.atomIdentifierMapModified(connectivity, atomMap);
        System.out.println(atomIdentifierMapModified + "atomIdentifierMapModified");
        ArrayList<ArrayList<Integer>> connectedAtoms = PDBReader.getConnectivity();
        System.out.println(connectedAtoms+ ": connectedAtom");
        ArrayList<ArrayList<Integer>> connectivityModified = PDBReader.getconnectivityModified(connectedAtoms);
        System.out.println(connectivityModified+ ": connectedAtomModified" );
        Map<Integer,String> atomMap = PDBReader.getAtomMap(connectedAtoms);
        System.out.println(atomMap + ": atomMap");
        HashMap<Integer, String> atomMapModified = PDBReader.getatomMapModified(atomMap);
        System.out.println(atomMapModified + ": atomMapModified");
        getUniqueAtoms(connectivityModified, atomMapModified);
        System.out.println(uniqueElements + "Set of Unique Elements");
        System.out.println(atomIdentifierMap +" atomIdentifierMap");
        System.out.println(atomIdentifierMapModified + " atomIdentifierMapModified");

        //List<int[]>valuePairs = firstLocationOfAtomPrinter(atomIdentifierMap);
        //System.out.println(Arrays.deepToString(valuePairs.toArray()) + "valuePairs");
        //List<int[]> listOfBonds = getBondList(connectivityModified);
        //List<int[]> listOfAngleModified = getAngleList(connectivityModified);
        //List<int[]> listOfTorsionModified = getTorsionList(connectivity);
        //System.out.println(Arrays.deepToString(listOfBonds.toArray())+ ": listOfBond");
        //Map<Integer, String> atomIdentifierMapMod = PDBReader.atomIdentifierMapModified(connectivityModified, atomMapModified);
        //List<int[]>dupletsSorted= bondSorter(listOfBonds, atomIdentifierMapMod);
       // List<int[]>tripletsSorted=angleSorter(listOfAngleModified, atomIdentifierMapMod);
      //  List<int[]>quadrupletsSorted=torsionSorter(listOfTorsionModified, atomIdentifierMapMod);
       // System.out.println(Arrays.deepToString(quadrupletsSorted.toArray()) + "quadrupletsSorted");
        //System.out.println( Arrays.deepToString(listOfTorsionModified.toArray()) + "torsionModified");
        //List<int[]>quadrupletsSorted=torsionSorter(listOfAngleModified, atomIdentifierMapMod);
       // Map<String[],List<int[]>> bondTypesMap= idenBondTypes(dupletsSorted, atomIdentifierMapMod);
        //Map<String[],List<int[]>> angleTypesMap= idenAngleTypes(tripletsSorted, atomIdentifierMapMod);
       // Map<String[],List<int[]>> torsionTypesMap= idenTorsionTypes(quadrupletsSorted, atomIdentifierMapMod);
        //System.out.println(atomIdentifierMapMod);
        //idenBondTypes(listOfBonds, atomIdentifierMapMod);
       // uniqueBondTypes(listOfBonds, atomIdentifierMap);

       // System.out.println(Arrays.deepToString(listOfAngleModified.toArray())+ ": listOfAngleModified");
        //idenAngleTypes(listOfAngleModified, atomIdentifierMapMod);
       // List<int[]> listOfTorsionModified = getTorsionList(connectivity);
        //System.out.println(Arrays.deepToString(listOfTorsionModified.toArray())+ " listOfTorsionModified");
        //Set<String> uniqueAtoms = uniqueElementAtomTypeConverter(connectivity, atomMap);
        //System.out.println(uniqueAtoms + " uniqueAtoms");
        //System.out.println(atomIdentifierMap.size());
        //List<List<AtomType>> result = new ArrayList<>();
       // System.out.println("Values");
        // System.out.println(molecule + "molecule");
        for(int i = 0; i < atomIdentifierMap.size(); i++) {
            String symbol = String.valueOf(atomIdentifierMap.get(i));
            int startIndex = symbol.indexOf("[") + 1;
            int endIndex = symbol.indexOf("]");
            String nameNew = symbol.substring(startIndex, endIndex);
            //atomTypeRequirementIndividuals.add(nameNew);
           // System.out.println(atomIdentifierMap.get(i));
            AtomType newName = AtomType.simple(nameNew);
            //System.out.println(nameNew+ " :name " +  " " + newName + " :atomtype");
            if (typeMapNew.containsKey(nameNew)) {
                typeNew = typeMapNew.get(nameNew);
               // System.out.println(typeNew + "typenew");
            } else {
                typeNew = newName;
                typeMapNew.put(nameNew, typeNew);
                //System.out.println(typeMapNew + " typemapnew");
            }
            //System.out.println(positions.get(i) + " positions");
            Vector position = positions.get(i+1);
            speciesBuilderNew.addAtom(typeNew, position,  "");
        }
       // System.out.println(typeMapNew + " typeMapNew");
        return speciesBuilderNew.build();
    }
   /* public static List<int[]> firstLocationOfAtomPrinter( Map<Integer, String> atomIdentifierMap){
        Map<String, Integer> firstLocations = new HashMap<>();
        for (Map.Entry<Integer, String> entry : atomIdentifierMap.entrySet()) {
            int location = entry.getKey();
            String atomType = String.valueOf(entry.getValue());

            if (!firstLocations.containsKey(atomType)) {
                firstLocations.put(atomType, location);
            }
        }
        System.out.println(firstLocations + "firstLocations");

        List<int[]> valuePairs = new ArrayList<>();
        for (Map.Entry<String, Integer> entry1 : firstLocations.entrySet()) {
            String key1 = entry1.getKey();
            Integer value1 = entry1.getValue();

            for (Map.Entry<String, Integer> entry2 : firstLocations.entrySet()) {
                String key2 = entry2.getKey();
                Integer value2 = entry2.getValue();
                if(value1 >= value2){
                    int[] pair = {value1, value2};
                    System.out.println("Pair: " + key1 + "=" + value1 + ", " + key2 + "=" + value2);
                    valuePairs.add(pair);
                }

            }
        }
      return valuePairs;
    }

    public static List<List<String>> getAllAtomPairs(Set<String> atomList) {
        List<List<String>> result = new ArrayList<>();

        // Convert the atomList set to an array to make indexing easier
        String[] atomArray = atomList.toArray(new String[0]);
        System.out.println(Arrays.toString(atomArray) + " atomArray");
        // Loop through each pair of atoms, skipping duplicates (i.e. (A,B) is the same as (B,A))
        for (int i = 0; i < atomArray.length; i++) {
            for (int j = i; j < atomArray.length; j++) {
                List<String> currentPair = new ArrayList<>();
                currentPair.add(atomArray[i]);
                currentPair.add(atomArray[j]);
                result.add(currentPair);
            }
        }
        return result;
    }
    public static  List<List<AtomType>> atomPairs(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        Set<String> uniqueElements = PDBReader.uniqueElementIdentifierWithInputs(connectivity, atomMap);

    }*/

    public static Set<String> uniqueElementIdentifierWithInputs(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
    atomIdentifier(connectivityModified, atomMapModified);
    Set<String> uniqueAtoms = new HashSet<>();
    for(int i =0; i<atomIdentifierMap.size(); i++){
        AtomType atomName = atomIdentifierMap.get(i);
        //System.out.println(uniqueAtoms + "UniqueAtomsWithInputs Loop1");
        uniqueAtoms.add(String.valueOf(atomName));
    }
        return uniqueAtoms;
    }
    public static Map<Integer,String> atomIdentifierMapModified (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified) {
        Map<Integer, AtomType> atomIdentifierMap = atomIdentifier(connectivityModified, atomMapModified);
        Map<Integer, String> modifiedAtomIdentifierMap = new HashMap<>();
        for (Map.Entry<Integer, AtomType> entry : atomIdentifierMap.entrySet()) {
            String value = entry.getValue().toString();
            value = value.replace("AtomType[", "").replace("]", "");
            modifiedAtomIdentifierMap.put(entry.getKey(), value);
        }
        //System.out.println(modifiedAtomIdentifierMap);
        return modifiedAtomIdentifierMap;
    }
    public static ArrayList<Integer> getBondList (ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> atomMap){
        return bondList;
    }

    public static ArrayList<Integer> setBondList(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        ArrayList<Integer> bondList = singleDoubleBondIdentifier(connectivityModified, atomMapModified);
        return bondList;
    }
    public static ArrayList<Integer> singleDoubleBondIdentifier(ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> atomMap) {
        System.out.println("Already visited");
        int valenceCarbon = 4, valenceOxygen = 2, valenceHydrogen = 1, valenceNitrogen =0, valenceSulfur = 0, valenceHalide =1;
        int bondRequired = 0;
        ArrayList<Integer> nitrogenBonds;
        ArrayList<Integer> sulfurBonds;
        System.out.println(atomMap);
        for (int i = 0; i < connectivity.size(); i++) {
            int atomArraySize = connectivity.get(i).size();
            //System.out.println(atomArraySize + " start atomarraysize" +  connectivity.get(i));
            String atomName = atomMap.get(connectivity.get(i).get(0));
            if (atomName.equals("H")) {
                if (valenceHydrogen == atomArraySize - 1) {
                    // System.out.println("The hydrogen atom " + i + " is satisfied with one bond");
                    bondList.add(1);
                } else {
                    bondRequired = valenceHydrogen - (atomArraySize - 1);
                    //System.out.println("The hydrogen atom " + i + " is UNSATISFIED with " + bondRequired + " bond required");
                    bondList.add(0);
                }
            }
            if (atomName.equals("CL") ||atomName.equals("BR") || atomName.equals("F")) {
                if (valenceHydrogen == atomArraySize - 1) {
                    // System.out.println("The hydrogen atom " + i + " is satisfied with one bond");
                    bondList.add(1);
                } else {
                    bondRequired = valenceHydrogen - (atomArraySize - 1);
                    //System.out.println("The hydrogen atom " + i + " is UNSATISFIED with " + bondRequired + " bond required");
                    bondList.add(0);
                }
            }

            if (atomName.equals("O")) {
                if (valenceOxygen == atomArraySize - 1) {
                    // System.out.println("The oxygen atom " + i + " is satisfied with one bond");
                    bondList.add(1);
                } else {
                    bondRequired = valenceOxygen - (atomArraySize - 1);
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
                // System.out.println(atomArraySize);
                if (valenceCarbon == atomArraySize - 1) {
                    //System.out.println("The carbon atom " + i + " is satisfied with one bond");
                    bondList.add(1);
                } else {
                    bondRequired = valenceCarbon - (atomArraySize - 1);
                    // System.out.println("The carbon atom " + i + " is UNSATISFIED with " + bondRequired + " bond required");
                    bondList.add(0);
                }
            }

            if (atomName.equals("N")){
                nitrogenBonds = connectivity.get(i);
                System.out.println(nitrogenBonds);
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
                    valenceNitrogen =4;
                    if( atomArraySize < 5){
                        bondList.add(0);
                    } else {
                        bondList.add(1);
                    }
                } else if (numOxy ==2) {
                    valenceNitrogen=4;
                    if( atomArraySize < 5){
                        bondList.add(0);
                    } else {
                        bondList.add(1);
                    }
                } else {
                    valenceNitrogen =3;
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
        System.out.println("Bonds are satisfied or not: "+ bondList);
        System.out.println(atomMap);
        ArrayList<Integer> connectivityArraylist = new ArrayList<>();
        ArrayList<Integer> connectivityArrayElementBondList = new ArrayList<>();
        int index;
        int e=0;
        while (e<=5){
            // System.out.println("Reached after while loop");
            for(int i =0; i<bondList.size(); i++){
                System.out.println(bondList.get(i) + " bondlist element " + i   );
                if (bondList.get(i) ==0){
                    System.out.println((i+1) + " " + bondList.get(i) + "i ");
                    System.out.println(connectivity.get(i) + " i" + i + "connectivity here");
                    System.out.println(connectivity.get(i) + " i" + i + "connectivity here");
                    connectivityArraylist = connectivity.get(i);
                    //System.out.println(connectivityArraylist + " connectivityarraylist");
                    for(int j =0; j<connectivityArraylist.size(); j++){
                        connectivityArrayElementBondList =connectivity.get(connectivityArraylist.get(j)-1);
                        index = connectivity.indexOf(connectivityArrayElementBondList);
                        if(bondList.get(index) == 0) {
                            //System.out.println(connectivityArrayElementBondList + " connectivityArrayElementBondList " + " " + bondList.get(index) + " bondlist element " + j + 1);
                            //System.out.println(connectivityArraylist + "Arraylist");
                            //System.out.println(connectivityArrayElementBondList + " arrayelmentlist");
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
                    // System.out.println(bondList + "bondlist");

                }
            }
            e++;
        }
        return bondList;
    }

    protected static  Map<Integer, AtomType> atomIdentifier(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        int counter =0;
        ArrayList <Integer> atomNumbers = new ArrayList<>();
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
                    IElement Carbon_3 = null;
                    AtomType C_3 = new AtomType(Carbon_3, "C_3");
                    //System.out.println("The atom " + (i) +" is C_3" );
                    atomIdentifierMap.put(i, C_3);
                }

                if (arrayListSize == 4) {
                    IElement Carbon_2 = null;
                    AtomType C_2 = new AtomType(Carbon_2, "C_2");
                    //System.out.println("The atom " + (i) +" is C_2 " );
                    atomIdentifierMap.put(i, C_2);
                }
                if (arrayListSize == 3) {
                    IElement Carbon_1 = null;
                    AtomType C_1 = new AtomType(Carbon_1, "C_1");
                    //System.out.println("The atom " + (i)+" is C_1 "  );
                    atomIdentifierMap.put(i, C_1);
                }
            } else if (retriveArrayFirstElementName.equals("CL")) {
                IElement Chlorine = null;
                AtomType Cl = new AtomType(Chlorine, "Cl");
                //System.out.println("The atom " + (i)+" is H " );
                atomIdentifierMap.put(i, Cl);
            } else if (retriveArrayFirstElementName.equals("I")) {
                IElement Iodine = null;
                AtomType I = new AtomType(Iodine, "I");
                //System.out.println("The atom " + (i)+" is H " );
                atomIdentifierMap.put(i, I);
            }else if (retriveArrayFirstElementName.equals("H")) {
                IElement Hydrogen = null;
                AtomType H = new AtomType(Hydrogen, "H");
                //System.out.println("The atom " + (i)+" is H " );
                atomIdentifierMap.put(i, H);
            }
            else if (retriveArrayFirstElementName.equals("O")) {
                int arrayListSize = connectivityModified.get(i).size();
                if(arrayListSize == 2){
                    IElement Oxygen_2 = null;
                    AtomType O_2 = new AtomType(Oxygen_2, "H");
                    //System.out.println("The atom " + (i)+" is O_2 " );
                    atomIdentifierMap.put(i, O_2);
                } else if (arrayListSize == 3) {
                    IElement Oxygen_3 = null;
                    AtomType O_3 = new AtomType(Oxygen_3, "H");
                    //System.out.println("The atom " + (i ) + " is O_3 ");
                    atomIdentifierMap.put(i, O_3);
                } else {
                    IElement Oxygen_Ar = null;
                    AtomType O_Ar = new AtomType(Oxygen_Ar, "O_Ar");
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
                    IElement Nitrogen_1 = null;
                    AtomType N_1 = new AtomType(Nitrogen_1, "N_1");
                    //System.out.println("The atom " + (i) +" is N_1 " );
                    atomIdentifierMap.put(i, N_1);

                } else if (arrayListSize ==4 && oxygenCounter>1) {
                    IElement Nitrogen_2 = null;
                    AtomType N_2 = new AtomType(Nitrogen_2, "N_2");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, N_2);
                } else {
                    IElement Nitrogen_3 = null;
                    AtomType N_3 = new AtomType(Nitrogen_3, "N_3");
                    //System.out.println("The atom " + (i) +" is N_3 " );
                    atomIdentifierMap.put(i, N_3);
                }

            } else if (retriveArrayFirstElementName.equals("S")) {
                //types are -S-, - -S=O =O, -S=O =O, -S=O =O -OH(-)
                int arrayListSize = connectivityModified.get(i).size();
                if(arrayListSize ==2){
                    IElement Sulfur_2= null;
                    AtomType S_2 = new AtomType(Sulfur_2, "S_2");
                    //System.out.println("The atom " + (i) +" is S_2 " );
                    atomIdentifierMap.put(i, S_2);
                }else {
                    IElement Sulfur_3 = null;
                    AtomType S_3 = new AtomType(Sulfur_3, "S_3");
                    //System.out.println("The atom " + (i) +" is S_3 " );
                    atomIdentifierMap.put(i, S_3);
                }

            } else if (retriveArrayFirstElementName.equals("P")){
                IElement Phosphorus_3 = null;
                AtomType P_3 = new AtomType(Phosphorus_3, "P_3");
                //System.out.println("The atom " + (i+1) +" is P_3 " );
                atomIdentifierMap.put(i, P_3);

            } else if (retriveArrayFirstElementName.equals("SI")) {
                IElement Silica = null;
                AtomType Si = new AtomType(Silica, "Si");
                // System.out.println("The atom " + (i+1) +" is Si " );
                atomIdentifierMap.put(i, Si);
                int arrayListSize = connectivityModified.get(i).size();
                 /*if(arrayListSize == 5){
                     for (int j = 1; j<5; j++){
                         String connectedElement = atomMap.get(connectivity.get(i).get(j));
                         System.out.println(connectedElement + "Send values");
                         atomIdentifierMap.put(i+1, AtomType.simple(connectedElement));
                     }
                 }*/

            }  else if (retriveArrayFirstElementName.equals("BR")) {
                IElement Bromine = null;
                AtomType Br = new AtomType(Bromine, "Br");
                //System.out.println("The atom " + (i) +" is Br " );
                atomIdentifierMap.put(i, Br);

            } else if (retriveArrayFirstElementName.equals("F")) {
                IElement Fluorine = null;
                AtomType F = new AtomType(Fluorine, "F");
                // System.out.println("The atom " + (i+1) +" is F " );
                atomIdentifierMap.put(i, F);

            } else {//Metal ions
                if (retriveArrayFirstElementName.equals("RH")) {
                    IElement Rhodium = null;
                    AtomType Rh = new AtomType(Rhodium, "Rh");
                    // System.out.println("The atom " + (i+1) +" is Rh " );
                    atomIdentifierMap.put(i, Rh);

                } else if (retriveArrayFirstElementName.equals("RU")) {
                    IElement Ruthenium = null;
                    AtomType Ru = new AtomType(Ruthenium, "Ru");
                    // System.out.println("The atom " + (i+1) +" is Ru " );
                    atomIdentifierMap.put(i, Ru);

                } else if (retriveArrayFirstElementName.equals("NI")) {
                    IElement Nickel = null;
                    AtomType Ni = new AtomType(Nickel, "Ni");
                    // System.out.println("The atom " + (i+1) +" is Ni " );
                    atomIdentifierMap.put(i, Ni);

                } else if (retriveArrayFirstElementName.equals("CU")) {
                    IElement Copper = null;
                    AtomType Cu = new AtomType(Copper, "Cu");
                    //System.out.println("The atom " + (i+1) +" is Cu " );
                    atomIdentifierMap.put(i, Cu);

                } else if (retriveArrayFirstElementName.equals("FE")) { // Not this
                    int arrayListSize = connectivityModified.get(i).size();
                    if (arrayListSize == 3) {
                        // connected to two elements
                        IElement Iron_2 = null;
                        AtomType Fe_2 = new AtomType(Iron_2, "Fe_2");
                        // System.out.println("The atom " + (i+1) +" is Fe_2 " );
                        atomIdentifierMap.put(i, Fe_2);

                    } else {
                        IElement Iron_3 = null;
                        AtomType Fe_3 = new AtomType(Iron_3, "Fe_3");
                        // System.out.println("The atom " + (i+1) +" is N_2 " );
                        atomIdentifierMap.put(i, Fe_3);

                    }
                    //System.out.println(atomValues);

                } else if (retriveArrayFirstElementName.equals("CO")) {
                    IElement Cobalt = null;
                    AtomType Co = new AtomType(Cobalt, "Co");
                    //System.out.println("The atom " + (i+1) +" is Co " );
                    atomIdentifierMap.put(i, Co);

                } else if (retriveArrayFirstElementName.equals("CR")) {
                    IElement Chromium = null;
                    AtomType Cr = new AtomType(Chromium, "Cr");
                    //System.out.println("The atom " + (i+1) +" is Cr " );
                    atomIdentifierMap.put(i, Cr);

                } else if (retriveArrayFirstElementName.equals("PD")) {
                    IElement Palladium = null;
                    AtomType Pd = new AtomType(Palladium, "Pd");
                    //System.out.println("The atom " + (i+1) +" is Pd " );
                    atomIdentifierMap.put(i, Pd);

                } else if (retriveArrayFirstElementName.equals("MO")) {
                    IElement Molybdenum = null;
                    AtomType Mo = new AtomType(Molybdenum, "Mo");
                    // System.out.println("The atom " + (i+1) +" is Mo " );
                    atomIdentifierMap.put(i, Mo);

                } else if (retriveArrayFirstElementName.equals("ZR")) {
                    IElement Zirconium = null;
                    AtomType Zr = new AtomType(Zirconium, "Zr");
                    // System.out.println("The atom " + (i+1) +" is Zr " );
                    atomIdentifierMap.put(i, Zr);

                } else if (retriveArrayFirstElementName.equals("V")) {
                    IElement Vanadium = null;
                    AtomType V = new AtomType(Vanadium, "V");
                    // System.out.println("The atom " + (i+1) +" is V " );
                    atomIdentifierMap.put(i, V);

                } else if (retriveArrayFirstElementName.equals("W")) {// Not this
                    IElement Tungsten = null;
                    AtomType W = new AtomType(Tungsten, "W");
                    //System.out.println("The atom " + (i+1) +" is W " );
                    atomIdentifierMap.put(i, W);

                } else if (retriveArrayFirstElementName.equals("MG")) {
                    IElement Magnesium = null;
                    AtomType Mg = new AtomType(Magnesium, "Mg");
                    // System.out.println("The atom " + (i+1) +" is Mg " );
                    atomIdentifierMap.put(i, Mg);

                } else if (retriveArrayFirstElementName.equals("ZN")) {
                    IElement Zinc = null;
                    AtomType Zn = new AtomType(Zinc, "Zn");
                    // System.out.println("The atom " + (i+1) +" is Zn " );
                    atomIdentifierMap.put(i, Zn);
                } else if (retriveArrayFirstElementName.equals("IR")) {
                    IElement Iridium = null;
                    AtomType Ir = new AtomType(Iridium, "Ir");
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
    public static Map<String, AtomType> getTypeMapNew(){
        return typeMapNew;
    }
    public static ArrayList<ArrayList<Integer>> getConnectivity(){
        //remove buildSpecies when running single file. When working along with
        // species= buildSpecies(confName);
        return connectivity;
    }

    public static ArrayList<ArrayList<Integer>> getConnectivityWithSpecies(String confName){
       buildSpecies(confName);
        return connectivity;
    }
    public static ArrayList<ArrayList<Integer>> getconnectivityModified (ArrayList<ArrayList<Integer>>connectivity){
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

    public static Map<Integer,String> getAtomMap( ArrayList<ArrayList<Integer>> connectivity){
        return atomMap;
    }
    public static  HashMap<Integer, String> getatomMapModified( Map<Integer, String> atomMap){
        for (Map.Entry<Integer, String> entry : atomMap.entrySet()) {
            atomMapModified.put(entry.getKey() - 1, entry.getValue());
        }
        return atomMapModified;
    }
    public static List<int[]> bondSorter(List<int[]> duplets,  Map<Integer, String> atomIdentifierMapModified){
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
        System.out.println(Arrays.deepToString(duplets.toArray()));
        System.out.println(Arrays.deepToString(dupletsActual.toArray()));
        System.out.println(atomIdentifierMapModified);
        //System.out.println(Arrays.deepToString(priorityList.toArray()));
        return dupletsActual;
    }
    public static List<int[]> angleSorter(List<int[]> triplets,  Map<Integer, String> atomIdentifierMapModified){
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
    public static List<int[]> torsionSorter(List<int[]> quadruplets,  Map<Integer, String> atomIdentifierMapModified){
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
    public static int priorityMapGenerator(String atomType){
        Map<String,Integer> priorityMap = new HashMap<>();
        priorityMap.put("Zr",1);
        priorityMap.put("C_3",2);
        priorityMap.put("C_2",3);
        priorityMap.put("C_1",4);
        priorityMap.put("H",5);
        priorityMap.put("Ph",6);
        priorityMap.put("N_3",7);
        priorityMap.put("N_2",8);
        priorityMap.put("N_1",9);
        priorityMap.put("O_3",10);
        priorityMap.put("O_2",11);
        priorityMap.put("O_1",12);
        priorityMap.put("I",13);
        priorityMap.put("Br",14);
        priorityMap.put("Cl",15);
        priorityMap.put("F",16);
        return priorityMap.get(atomType);
    }
    public static Map<String[],List<int[]>> idenBondTypes(List<int[]> duplets, Map<Integer, String> atomIdentifierMapModified){
        System.out.println(Arrays.deepToString(duplets.toArray()));
        Map<String[],List<int[]>> bondTypesMap = new HashMap<>();
        System.out.println(atomIdentifierMapModified);
        ArrayList<String[]> bondTypes = new ArrayList<>();
        System.out.println("\nThe printing of the bond types");
        //System.out.println("The printing of the bond types");
        //form bonded pairs of atoms
        for(int i =0; i< duplets.size(); i++){
            int firstElement = duplets.get(i)[0];
            int secondElement = duplets.get(i)[1];
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
                            // System.out.println("not equal " + newBond[j]);
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
                bondTypes.add(newBond);
            }
        }
        System.out.println(Arrays.deepToString(bondTypes.toArray()) + " bondTypes");
        //identify the pairs of individual types
        for(int i=0; i<bondTypes.size();i++){
            //System.out.println(bondTypes.get(i)[0] + " " + bondTypes.get(i)[1]);
            List<int[]> dupletsNew = new ArrayList<>();
            for(int k =0; k< duplets.size();k++){
                int firstElement = duplets.get(k)[0];
                int secondElement = duplets.get(k)[1];
                String atomOne = String.valueOf(atomIdentifierMapModified.get(firstElement));
                String atomTwo = String.valueOf(atomIdentifierMapModified.get(secondElement));
                // System.out.println(firstElement + " " + atomOneActual +" "+ secondElement + " " + atomTwoActual);
                if((bondTypes.get(i)[0].equals(atomOne) && bondTypes.get(i)[1].equals(atomTwo))||(bondTypes.get(i)[0].equals(atomTwo) && bondTypes.get(i)[1].equals(atomOne)) ){
                    int[] sameBond = {firstElement, secondElement};
                    dupletsNew.add(sameBond);
                }
                //System.out.println(Arrays.deepToString(dupletsNew.toArray()) + " dupletsNew After");
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
        System.out.println("\n");
        return bondTypesMap;
    }

    public static Map<String[],List<int[]>> idenAngleTypes( List<int[]> triplets,  Map<Integer, String> atomIdentifierMapModified){
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
    public static Map<String[],List<int[]>> idenTorsionTypes( List<int[]> quadruplets,  Map<Integer, String> atomIdentifierMapModified){
         //System.out.println(Arrays.deepToString(quadruplets.toArray()) + "quadruplets");
        Map<String[],List<int[]>> torsionTypesMap = new HashMap<>();
        //System.out.println(atomIdentifierMapModified);
        ArrayList<String[]> torsionTypes = new ArrayList<>();
        System.out.println("The printing of the torsion types");
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
        System.out.println(Arrays.deepToString(torsionTypes.toArray()) + " torsionTypes");
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


    public static List<int[]> getBondedAtomList (ArrayList<ArrayList<Integer>> connectivityModified){
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
        List<int[]> intListArray = new ArrayList<>();
        for (List<Integer> innerList : listOfBonds) {
            int[] array = innerList.stream().mapToInt(Integer::intValue).toArray();
            intListArray.add(array);
        }
        return intListArray;
    }
    public static ArrayList<ArrayList<ArrayList<Integer>>> generateCombinationsForNestedList(ArrayList<ArrayList<Integer>> inputList) {
        ArrayList<ArrayList<ArrayList<Integer>>> outputList = new ArrayList<>();
        for (ArrayList<Integer> innerList : inputList) {
            ArrayList<ArrayList<Integer>> innerOutputList = generateCombinations(innerList);
            outputList.add(innerOutputList);
        }
        return outputList;
    }
    public static ArrayList<ArrayList<Integer>> generateCombinations(List<Integer> inputList) {
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

    public static List<int[]> getAngleList(ArrayList<ArrayList<Integer>> connectivityModified) {
        ArrayList<ArrayList<ArrayList<Integer>>> arrayListOfAngles = generateCombinationsForNestedList(connectivityModified);
        List<int[]> result = new ArrayList<>();

        for (ArrayList<ArrayList<Integer>> subList1 : arrayListOfAngles) {
            for (ArrayList<Integer> subList2 : subList1) {
                int[] arr = subList2.stream().mapToInt(Integer::intValue).toArray();
                result.add(arr);
            }
        }
        return result;
    }
    public static ArrayList<Integer> getlistOfTorsions(ArrayList<ArrayList<Integer>> connectivityModified){
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
    public static List<List<Integer>> setTorsionPairs(ArrayList<ArrayList<Integer>> connectivityModified){
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

    public static List<int[]> getTorsionList(ArrayList<ArrayList<Integer>> connectivity) {
        List<List<Integer>> listOfTorsion = setTorsionPairs(connectivity);
        List<int[]> listOfTorsionModified = listOfTorsion.stream()
                .map(sublist -> new int[]{sublist.get(0) - 1, sublist.get(1) - 1, sublist.get(2) - 1, sublist.get(3) - 1})
                .collect(Collectors.toList());
        //System.out.println(listOfTorsionModified);
        return listOfTorsionModified;

    }


    public static Set<String> uniqueElementIdentifier(){
        atomIdentifier(connectivityModified, atomMapModified);
        //System.out.println(atomIdentifierMap);
        Set<String> uniqueAtoms = new HashSet<>();
        //System.out.println(atomIdentifierMap + "in UniqueElement Identifier @761");
       /*for(int i =0; i< atomIdentifierMap.size(); i++){
            System.out.println(atomIdentifierMap.get(i));
            //AtomType atomName = atomIdentifierMap.get(i);
            //System.out.println(atomName);
        }*/
        for(int i =0; i<atomIdentifierMap.size(); i++){
            AtomType atomName = atomIdentifierMap.get(i);
            //System.out.println(atomName + "uniqueElement Identifier");
            uniqueAtoms.add(String.valueOf(atomName));
        }
        // System.out.println(uniqueAtoms + "UniqueAtoms");
        return uniqueAtoms;
    }

    public static Set<String> getUniqueAtoms(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        Set<String> uniqueAtomsWithInputs= uniqueElementIdentifier();
        //System.out.println(uniqueAtomsWithInputs + "uniqueAtom WithInputs");
        Set<String> result = new HashSet<>();
        for (String atomTypeString : uniqueAtomsWithInputs) {
            String[] parts = atomTypeString.split("\\[|\\]");
            //System.out.println(Arrays.toString(parts));
            result.add(parts[1]);
        }
        System.out.println(result + "results");
        return result;
    }

    public static double getMolecularWeight(){
        //String  molecularWeightString = String.format("%.4f", molecularWeight);
        //System.out.println( "Molecular Weight: "+molecularWeightString);
        return molecularWeight;
    }

    public static void getUniqueAtomsForOtherClasses(String confName){
        ArrayList<ArrayList<Integer>>connectivity= getConnectivityWithSpecies(confName);
        ArrayList<ArrayList<Integer>>connectivityModified = getconnectivityModified(connectivity);
        Map<Integer,String> atomMap = getAtomMap(connectivity);
        HashMap<Integer, String> atomMapModified = getatomMapModified(atomMap);
        //System.out.println(connectivityModified);
        //System.out.println(atomMapModified);
    }


    public static Map<String, double[]> atomicPotMap(){
        Map<String, double[]> atomicPot = new HashMap<>();
        Set<String> uniqueElement= uniqueElementIdentifier();
        Iterator<String> iterator = uniqueElement.iterator();
        while (iterator.hasNext()){
            String str = iterator.next();
            //System.out.println(str + "str");
            String elementName = str.substring(str.indexOf("[") + 1, str.indexOf("]"));
            //System.out.println( str+ " " +elementName);
            double[] atomicPotValues = atomicPot(elementName);
            //System.out.println(str + " " + Arrays.toString(atomicPotValues));
            //atomicPotMap.put(str, atomicPotValues);
            atomicPot.put(elementName, atomicPotValues);
            System.out.println(elementName + " " + atomicPotValues);
        }
        return atomicPot;
    }

    public static int doubleAtomIdentifierForTorsion(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        int numCarbon = 0;
        for(int i =0; i<connectivityModified.size(); i++){
            String atomName = atomMapModified.get(connectivityModified.get(i).get(0));{
                if(atomName.equals("C")){
                    System.out.println(connectivityModified.get(i) + "connectivityModified" + " "+ atomName + " inside");
                    numCarbon++;
                }
            }
        }
        int y = 0;
        if(numCarbon % 2 ==0){
            y = 1;
        } else {
            y=0;
        }
        return y;
    }

    public static double [] atomicPot (String atomtype){
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
        atomicConstant.put("N_2", new double[]{0.699, 120.0, 3.66, 0.069, 13.407, 2.544,6.899,2.0});
        atomicConstant.put("N_3", new double[]{0.7, 106.7, 3.66, 0.069, 13.407, 2.544, 6.899,0.45});
        atomicConstant.put("N_1", new double[]{0.685, 111.2, 3.66, 0.069, 13.407, 2.544,6.899,2.0});
        atomicConstant.put("O_2", new double[]{0.68, 110.0, 3.5, 0.06, 14.085, 2.3, 8.741, 2.0});
        atomicConstant.put("O_3", new double[]{0.658, 104.51, 3.5, 0.06, 14.085, 2.3,8.74,0.018});
        atomicConstant.put("O_1", new double[]{0.639, 180.0, 3.5, 0.06, 14.085, 2.3, 8.741,2.0});
        atomicConstant.put("S_2", new double[]{1.077, 92.2, 4.035, 0.274, 13.969, 2.703,6.928,1.25});
        atomicConstant.put("S_3", new double[]{1.064, 92.1, 4.035, 0.274, 13.969, 2.703,6.928,0.484});
        atomicConstant.put("S+4", new double[]{1.049, 103.2, 4.035, 0.274, 13.969, 2.703,6.928,0.484});
        atomicConstant.put("S+6", new double[]{1.027, 109.47, 4.035, 0.274, 13.969, 2.703,6.928,0.484});
        atomicConstant.put("P_3", new double[]{1.101, 93.8, 4.147, 0.305, 13.072, 2.863,5.463,2.4});
        //atomicConstant.put("P+5", new double[]{1.056, 109.47, 4.147, 0.305, 13.072, 2.863,5.463});
        atomicConstant.put("H", new double[]{0.354, 180.0, 2.886, 0.044, 12.0, 0.71,4.528, 0.0});
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
        double [] sample = atomicConstant.get(atomtype);
        return sample;
    }
    public static double getMolecularWeight (String symbol){
        // Define the atomic weights of the elements
        switch (symbol){
            case "H":
                return 1.0079;
            case "C":
                return 12.011;
            case "N":
                return 14.0067;
            case "O":
                return 15.9994;
            case "F":
                return 18.998;
            case "Na":
                return 22.9897;
            case "Mg":
                return 24.305;
            case "Al":
                return 26.9815;
            case "Si":
                return 28.0855;
            case "P":
                return 30.9738;
            case "S":
                return 32.06;
            case "Cl":
                return 35.453;
            case "K":
                return 	39.0983;
            case "Ti":
                return 47.90;
            case "V":
                return 50.9415;
            case "Fe":
                return 	55.847;
            case "Cr":
                return 51.996;
            case "Mn":
                return 54.9380;
            case "Ni":
                return 	58.70;
            case "Co":
                return 58.9332;
            case "Cu":
                return 63.546;
            case "Zn":
                return 65.38;
            case "Br":
                return 79.904;
            case "Zr":
                return 91.22;
            case "Mo":
                return 95.94;
            case "Ru":
                return 101.07;
            case "Rh":
                return 102.9055;
            case "Pd":
                return 	106.4;
            case "Sb":
                return 121.75;
            case "Re":
                return 186.207;
            case "Ir":
                return 192.22;
            case "Sn":
                return 118.69;
            case "I":
                return 126.9045;
            default:
                return 0.00;
        }
    }

    public static double [] electronicParam (String atomtype){
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
    public static void structureCompare(String struc1, String struc2){
        String confName = struc1;
        String conf2Name = struc2;
        species = getSpeciesNew(confName);
        ArrayList<ArrayList<Integer>> benzeneConnectivity = connectivity;
        Map<Integer, String> benzeneAtomMap = atomMap;
        Map<Integer, AtomType> benzeneAtomIdentifier = atomIdentifierMap;
        System.out.println(benzeneAtomMap+ " benzene");
        System.out.println(benzeneAtomIdentifier+ " benzene");
        System.out.println(species + "species NEW");
        bondList.clear();
        connectivity.clear();
        atomMap.clear();
        atomIdentifierMap.clear();
        connectivityModified.clear();

        //bondList = null;
        //connectivity = null;
       // atomMap = null;
        speciesNew = getSpeciesNew(conf2Name);
        ArrayList<ArrayList<Integer>> chlorobenzeneConnectivity = connectivity;
        Map<Integer, String> chlorobenzeneAtomMap = atomMap;
        Map<Integer, AtomType> chlorAtomIdentifier = atomIdentifierMap;
        System.out.println(chlorobenzeneAtomMap + " chlorbenzene");
        System.out.println(chlorAtomIdentifier + " chlorbenzene");
        System.out.println(speciesNew + "species NEW");
        System.out.println(species + " species");
        boolean isBenzeneInChlorobenzene = chlorobenzeneConnectivity.containsAll(benzeneConnectivity)
                && benzeneAtomIdentifier.entrySet().containsAll(chlorAtomIdentifier.entrySet());
        File file = new File(conf2Name);
        String confFileName = file.getName();
        if (isBenzeneInChlorobenzene) {
            System.out.println("Benzene structure is present in " + confFileName);
        } else {
            System.out.println("Benzene structure is not present in " + confFileName);
        }
    }

    public static ArrayList<ArrayList<Integer>> getConnectivityWithoutRunning(){
        return connectivity;
    }

    public static ArrayList<ArrayList<Integer>> getConnectivityModifiedWithoutRunning(){
        return connectivityModified;
    }
    public static Map<Integer,String> getAtomMapWithoutRunning(){
        return atomMap;
    }

    public static HashMap<Integer, String> getAtomMapModifiedWithoutRunning(){
        return atomMapModified;
    }



    public static void main(String[] args) {
        //structureCompare("F:/benzene", "F:/ethane");
        /*String confName = "F:/ethane";
        speciesNew =getSpeciesNew(confName);
        System.out.println(speciesNew + "species NEW");
        String confNameNew = "F:/benzene";
        species =getSpeciesNew(confNameNew);
        System.out.println(species + "species NEW");*/

    }
}
