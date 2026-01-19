package etomica.potential.GAFF;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.chem.elements.*;
import etomica.molecule.IMolecule;
import etomica.potential.UFF.PDBReader;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class PDBReaderGAFF {
    ISpecies species;
    Map<Integer, Integer> assignCharges = new HashMap<>();
    Map<Integer, AtomType> atomIdentifierGAFF = new HashMap<>();
    static Map<Integer, String> modifiedAtomIdentifierMap = new HashMap<>();
    public static Map<String, AtomType> elementReceiverMap = new HashMap<>();
    public static Map<Integer,Vector> positions = new HashMap<>();
    public Map<Double, Double> chargeMap = new HashMap<>();

    public Map<Double, Double> readChargeFile(String confName) {
        String fileName = confName+".pdb";
        FileReader fileReader;

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
                line = line.trim();
                if (line.length() < 6) {
                    throw new RuntimeException("Incorrect") ;
                }
                if (line.substring(0, 6).equals("HETATM") || line.substring(0,4).equals("ATOM")) {
                    //coordinates of atom and create atom
                    double x = Double.parseDouble(line.substring(30, 38));
                    double y = Double.parseDouble(line.substring(38, 46));
                    chargeMap.put(x, y);
                }
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        return chargeMap;
    }

    public ISpecies getSpecies (String confName, String chargeName){
        double massSum = 0;
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        SpeciesBuilder speciesBuilderNewMod =  new SpeciesBuilder(Space3D.getInstance());
        AtomType typeNew;
        PDBReader pdbReader = new PDBReader();
        pdbReader.readPDBFile(confName);
        readChargeFile(chargeName);
        ArrayList<ArrayList<Integer>> connectivity = pdbReader.getConnectivity();
        Map<Integer, String> atomMap = pdbReader.getAtomMap();
        ArrayList<Integer> bondList =pdbReader.setBondList(connectivity, atomMap);
        //System.out.println(bondList + " bondList");
        System.out.println(atomMap+ "  " + atomMap.size());
        Map<Integer, String> atomIdentifierMapMod = atomIdentifierMapModifiedGAFF(connectivity, atomMap);
        // System.out.println(atomIdentifierMapModified + "atomIdentifierMapModified");
        ArrayList<ArrayList<Integer>> connectedAtoms =pdbReader.getConnectivity();
        //System.out.println(connectedAtoms+ ": connectedAtom");
        ArrayList<ArrayList<Integer>> connectivityModified = pdbReader.getconnectivityModified(connectedAtoms);
        // System.out.println(connectivityModified+ ": connectedAtomModified" );
        // System.out.println(atomMap + ": atomMap");
        HashMap<Integer, String> atomMapModified = pdbReader.getatomMapModified(atomMap);
        // System.out.println(atomMapModified + ": atomMapModified");
        positions = pdbReader.getPositions();
        System.out.println(connectedAtoms +" O");
        System.out.println(connectivityModified+" O");
        System.out.println(atomMapModified+" O");
        System.out.println(atomIdentifierMapMod+" O");
        ArrayList<Integer> bondsNum = pdbReader.bondsAmongAtoms();
        System.out.println(bondsNum + " bond Amongatoms");
        List<int[]> duplets =  pdbReader.getBondedAtomList(connectivityModified);
        List<int[]> listOfAngleModified =  pdbReader.getAngleList(connectivityModified);
        List<int[]> dupletsSorted = pdbReader.bondSorter(duplets, atomIdentifierMapMod);
        //System.out.println(Arrays.deepToString(dupletsSorted.toArray()) + " duplets Sorted");
        //System.out.println(Arrays.deepToString(duplets.toArray())+ ": listOfBonds");
        List<int[]> tripletsSorted = pdbReader.angleSorter(listOfAngleModified, atomIdentifierMapMod);
        List<int[]> listOfTorsionModified =  pdbReader.getTorsionList(connectivity);
        //System.out.println( Arrays.deepToString(listOfTorsionModified.toArray()) + "torsionModified");
        List<int[]> quadrupletsSorted = pdbReader.torsionSorter(listOfTorsionModified, atomIdentifierMapMod);
        // System.out.println(Arrays.deepToString(quadrupletsSorted.toArray()) + "quadrupletsSorted");
        Map<String[],List<int[]>> torsionTypesMap=  pdbReader.idenTorsionTypes(quadrupletsSorted, atomIdentifierMapMod);
        List<int[]> listOfInversions = pdbReader.idenInversions(tripletsSorted);
        //System.out.println(Arrays.deepToString(listOfInversions.toArray()) + " in Main" );
        System.out.println(pdbReader.positions + " Here positions");
        Space space = Space3D.getInstance();
        Vector center = space.makeVector();
        Vector dr = Vector.d(center.getD());
        System.out.println(atomIdentifierMapMod);
        Map<String, AtomType> typeMapNew = pdbReader.getTypeMapNew();
        for(int i = 0; i < atomIdentifierMapMod.size(); i++) {
            String symbol = String.valueOf(atomIdentifierMapMod.get(i));
            int startIndex = symbol.indexOf("[") + 1;
            int endIndex = symbol.indexOf("]");
            String nameNew = symbol.substring(startIndex, endIndex);
            // System.out.println(atomIdentifierMap.get(i));
            AtomType newName = returnElementGAFF(nameNew);
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
        for(int i=0; i<atomIdentifierMapMod.size(); i++){
            IAtom a = children.get(i);
            String symbol = String.valueOf(atomIdentifierMapMod.get(i));
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

    protected Map<Integer, AtomType> atomIdentifierGAFF(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        for(int i =0; i<connectivityModified.size(); i++){
            Integer retriveArrayFirstElementNumber = connectivityModified.get(i).get(0);
            String retriveArrayFirstElementName = atomMapModified.get(retriveArrayFirstElementNumber);
            int arrayListSize = connectivityModified.get(i).size();
            int numCarbon=0, numOxy=0, numNitro=0, numHydro=0;
            ArrayList<Integer> innerList=connectivityModified.get(i);
            ArrayList<Integer> innerAtoms= numAtomCounter(connectivityModified, atomMapModified, innerList);
            numOxy=innerAtoms.get(2);
            if(retriveArrayFirstElementName.equals("C")){
                if(arrayListSize == 5){
                    AtomType CT = new AtomType(Carbon.INSTANCE, "CT");
                    atomIdentifierGAFF.put(i, CT);
                } else if (arrayListSize == 4) {
                    if(numOxy == 0){
                        // General alkene carbon
                        AtomType CM = new AtomType(Carbon.INSTANCE, "CM");
                        atomIdentifierGAFF.put(i, CM);
                    } else  {
                        //Acid carbon
                        AtomType C = new AtomType(Carbon.INSTANCE, "C");
                        atomIdentifierGAFF.put(i, C);
                    }
                }
            } else if (retriveArrayFirstElementName.equals("H")) {
                Integer atomConnected = connectivityModified.get(i).get(1);
                String atomNameConnected = atomMapModified.get(atomConnected);
                if(atomNameConnected.equals("C")){
                    ArrayList<Integer> atomSecondaryConnected = connectivityModified.get(atomConnected);
                    ArrayList<Integer> atomNameSecondaryConnected = secondaryAtomCounter(connectivityModified, atomMapModified, atomSecondaryConnected);
                    if(atomNameSecondaryConnected.get(2) == 0 && atomNameSecondaryConnected.get(3) == 0){
                        //Connected to Carbon with one electronegative group
                        AtomType HC = new AtomType(Hydrogen.INSTANCE, "HC");
                        atomIdentifierGAFF.put(i, HC);
                    } else if(atomNameSecondaryConnected.get(2) == 1 || atomNameSecondaryConnected.get(3) == 1){
                        //Connected to Carbon with one electronegative group
                        AtomType H1 = new AtomType(Hydrogen.INSTANCE, "H1");
                        atomIdentifierGAFF.put(i, H1);
                    } else if (atomNameSecondaryConnected.get(2) == 2 || atomNameSecondaryConnected.get(3) == 2) {
                        AtomType H2 = new AtomType(Hydrogen.INSTANCE, "H2");
                        atomIdentifierGAFF.put(i, H2);
                    }
                } else if (atomNameConnected.equals("N")) {
                    AtomType H = new AtomType(Hydrogen.INSTANCE, "H");
                    atomIdentifierGAFF.put(i, H);
                } else if (atomNameConnected.equals("O")) {
                    AtomType HO = new AtomType(Hydrogen.INSTANCE, "HO");
                    atomIdentifierGAFF.put(i, HO);
                }
            }else if (retriveArrayFirstElementName.equals("N")) {
                if(arrayListSize ==4){
                    AtomType N3 = new AtomType(Nitrogen.INSTANCE, "N3");
                    atomIdentifierGAFF.put(i, N3);
                }
            }else if (retriveArrayFirstElementName.equals("O")) {
                Integer atomConnected = connectivityModified.get(i).get(1);
                String atomNameConnected = atomMapModified.get(atomConnected);
                if(atomNameConnected.equals("H")){
                    //alcohol
                    AtomType OH = new AtomType(Oxygen.INSTANCE, "OH");
                    atomIdentifierGAFF.put(i, OH);
                } else if (innerAtoms.get(0) == 2) {
                    //ester
                    AtomType OS = new AtomType(Oxygen.INSTANCE, "OS");
                    atomIdentifierGAFF.put(i, OS);
                } else if (arrayListSize ==2) {
                    AtomType OH = new AtomType(Oxygen.INSTANCE, "OH");
                    atomIdentifierGAFF.put(i, OH);
                }else {
                    ArrayList<Integer> atomSecondaryConnected = connectivityModified.get(atomConnected);
                    ArrayList<Integer> atomNameSecondaryConnected = secondaryAtomCounter(connectivityModified, atomMapModified, atomSecondaryConnected);
                    if(atomNameSecondaryConnected.get(3) == 1 && atomNameSecondaryConnected.get(2) ==1){
                        AtomType O = new AtomType(Oxygen.INSTANCE, "O");
                        atomIdentifierGAFF.put(i, O);
                    }
                }

            }else {
                throw new RuntimeException("Atom doesnot exist");
            }
        }
        return atomIdentifierGAFF;
    }
    public static AtomType returnElementGAFF(String elementName){
        elementReceiverMap.put("CT", new AtomType(Carbon.INSTANCE,"CT"));
        elementReceiverMap.put("CM", new AtomType(Carbon.INSTANCE,"CM"));
        elementReceiverMap.put("C", new AtomType(Carbon.INSTANCE,"C"));
        elementReceiverMap.put("H1", new AtomType(Carbon.INSTANCE,"H1"));
        elementReceiverMap.put("HC", new AtomType(Carbon.INSTANCE,"HC"));
        elementReceiverMap.put("H2", new AtomType(Carbon.INSTANCE,"H2"));
        elementReceiverMap.put("H", new AtomType(Carbon.INSTANCE,"H"));
        elementReceiverMap.put("HO", new AtomType(Carbon.INSTANCE,"HO"));
        elementReceiverMap.put("N3", new AtomType(Carbon.INSTANCE,"N3"));
        elementReceiverMap.put("OH", new AtomType(Carbon.INSTANCE,"OH"));
        elementReceiverMap.put("OS", new AtomType(Carbon.INSTANCE,"OS"));
        elementReceiverMap.put("O", new AtomType(Carbon.INSTANCE,"O"));
        return elementReceiverMap.get(elementName);
    }

    public static Map<Integer,String> atomIdentifierMapModifiedGAFF (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        PDBReaderGAFF readerGAFF = new PDBReaderGAFF();
        Map<Integer, AtomType> atomIdentifierMap = readerGAFF.atomIdentifierGAFF(connectivityModified, atomMapModified);
        for (Map.Entry<Integer, AtomType> entry : atomIdentifierMap.entrySet()) {
            String value = entry.getValue().toString();

            value = value.replace("AtomType[", "").replace("]", "");
            modifiedAtomIdentifierMap.put(entry.getKey(), value);
        }
        System.out.println(modifiedAtomIdentifierMap + " modified Identifier");
        return modifiedAtomIdentifierMap;
    }

    public ArrayList<Integer> numAtomCounter (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified, ArrayList<Integer> innerList){
        int numCarbon=0, numOxy=0, numNitro=0, numHydro=0;
        ArrayList<Integer> numAtomCounter = new ArrayList<>();
        for(int j =0; j< innerList.size(); j++){
            String innerAtomName = atomMapModified.get(innerList.get(j));
            if(innerAtomName.equals("O")){
                numOxy++;
            } else if (innerAtomName.equals("N")) {
                numNitro++;
            } else if (innerAtomName.equals("C")) {
                numCarbon++;
            } else if (innerAtomName.equals("H")) {
                numHydro++;
            } else {
                throw new RuntimeException("Atom does not exist");
            }
        }
        numAtomCounter.add(0,numCarbon);
        numAtomCounter.add(1,numHydro);
        numAtomCounter.add(2,numOxy);
        numAtomCounter.add(3,numNitro);
        return numAtomCounter;
    }
    public ArrayList<Integer> secondaryAtomCounter (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified, ArrayList<Integer> innerList){
        return numAtomCounter(connectivityModified, atomMapModified, innerList);
    }
    class InvalidAtomException extends Exception {
        public InvalidAtomException(String message) {
            super(message);
        }
    }

    public static void main(String[] args) throws InvalidAtomException {
        String confName ="";
        String chargeName = "";
        PDBReaderGAFF readerGAFF = new PDBReaderGAFF();
        ISpecies species = readerGAFF.getSpecies(confName, chargeName);

    }
}
