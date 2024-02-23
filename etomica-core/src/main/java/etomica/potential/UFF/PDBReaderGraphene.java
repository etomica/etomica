package etomica.potential.UFF;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Nitrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

//This class provide bonds within graphene
public class PDBReaderGraphene {
    static ISpecies species;
    double massSum;
    public Map<Integer, String> atomMapGraphene = new HashMap<>();
    public Map<Integer, Vector> positions = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectivityGrapehene = new ArrayList<>();
    public Map<Integer, String> atomMap = new HashMap<>();
    public Map<Integer, ArrayList<Integer>> connectivityMap = new HashMap<>();
    public  ArrayList<Integer> missingElement = new ArrayList<>();
    public Map<Integer, String> atomMapMainStruc = new HashMap<>();
    public Map<Integer, String> modifiedAtomIdentifierMapGraphene = new HashMap<>();
    public Map<Integer, Vector> positionsGraphene = new HashMap<>();
    public Map<String, AtomType> elementReceiverMap = new HashMap<>();
    public Map<String, AtomType> typeMap = new HashMap<>();
    public HashSet<String> uniqueElements = new HashSet<>();

    public void readPDBGrapheneFile(String confName) {
        String fileName = confName + ".pdb";
        FileReader fileReader;
        ArrayList<Integer> currentAtomList = null;
        try {
            fileReader = new FileReader(fileName);
        } catch (IOException e) {
            throw new RuntimeException("Cannot open " + fileName + ", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            String line;
            while ((line = bufReader.readLine()) != null) {
                parseLineReader(line, atomMap,atomMapMainStruc, positions, uniqueElements);
            }
            fileReader.close();
        } catch (IOException e) {
            throw new RuntimeException("Problem reading from " + fileName + ", caught IOException: " + e.getMessage());
        }

        //String  molecularWeightString = String.format("%.4f", molecularWeight);
        // System.out.println( "Molecular Weight: "+molecularWeightString);
        //System.out.println(typeMap + "typeMap");
       // System.out.println(uniqueElements);
        //System.out.println(atomMap + " atomMap");
        System.out.println("Building done");
        //System.out.println(positions);
        //System.out.println("\n " + atomMapMainStruc);
        Map<Integer, String> sortedMapMainStruc = new TreeMap<>(atomMapMainStruc);
       // System.out.println("\n "+sortedMapMainStruc);
        Map<Integer, String> sortedMap = new TreeMap<>(atomMap);
        //System.out.println("\n "+sortedMap);
        int expectedKey = 1;
        for (int i =0; i<sortedMap.size(); i++){
            if(sortedMap.get(i )== null){
                missingElement.add(i);
            }
        }
        for (int key : sortedMap.keySet()) {
            while (expectedKey < key) {
                missingElement.add(expectedKey);
                throw new IllegalStateException("Integer " + expectedKey + " is missing in the sorted map.");
            }
            expectedKey++;
        }
        System.out.println(missingElement);
        Vector dr = new Vector3D();
        for(int i=0; i<positions.size(); i++){
            ArrayList<Integer> newInt = new ArrayList<>();
            connectivityMap.put(i, newInt);
        }
        for (int i=0; i<positions.size(); i++){
            Vector newPos = new Vector3D();
            Vector iPosition = positions.get(i);
            for (int j=0; j<positions.size(); j++){
                ArrayList<Integer>  newInt = new ArrayList<>();
                ArrayList<Integer> newInt2 = new ArrayList<>();
                Vector jPosition = positions.get(j);
                newPos.E(iPosition);
                if(i<j){
                    //System.out.println(i + " "+ j);
                    double value = newPos.Mv1Squared(jPosition);
                    if(2.0101<value && value<2.0109){
                        newInt = connectivityMap.get(i);
                        newInt.add(j);
                        connectivityMap.put(i, newInt);
                        newInt2 = connectivityMap.get(j);
                        newInt2.add(i);
                        connectivityMap.put(j, newInt2);
                    }
                }
            }
        }
        System.out.println(connectivityMap);
        //System.out.println(positions);
        System.out.println(missingElement);;
    }
    protected void parseLineReader(String line, Map<Integer, String> atomMap,Map<Integer, String> atomMapMainStruc, Map<Integer, Vector> positions, HashSet<String> uniqueElements) {
        line = line.trim();
        if (line.length() < 6) {
            return;
        }
        if (line.substring(0, 6).equals("HETATM") || line.substring(0,4).equals("ATOM")) {
            //coordinates of atom and create atom
            double x = Double.parseDouble(line.substring(30, 38));
            double y = Double.parseDouble(line.substring(38, 46));
            double z = Double.parseDouble(line.substring(46, 54));
            String atom = String.valueOf(line.substring(13,15)).trim();
            Vector positn = Vector.of(x, y, z);
            String symbol = line.substring(12, 16).trim();
            int atomNumber = Integer.parseInt(line.substring(8,11).trim());
            positions.put(atomNumber, positn);
            uniqueElements.add(atom);
            if (atomMap.containsKey(atomNumber)){
                //System.out.println("Error");
            } else {
                // tempMap.put(atomNumber, typeElement);
                atomMap.put(atomNumber, symbol);
                if(z==0){
                    atomMapMainStruc.put(atomNumber, symbol);
                }
            }
        }
    }


    public Map<Integer, String> getModifiedAtomIdentifierMap(Map<Integer, String> atomMapModified) {
        HashSet<String> uniquesAtoms = new HashSet<>();
        for (String str : atomMapModified.values()) {
            uniquesAtoms.add(str);
        }
        System.out.println(uniquesAtoms);
        int n = 0;
        if (atomMapModified.get(0) == null) {
            n = 1;
        }
        for (String str : atomMapModified.values()) {
            n++;
            if (str == null) {
                throw new RuntimeException("Does not exist " + n);
            }
        }
        Map<Integer, AtomType> atomIdentifierMapGAFF = atomIdentifierGraphene(atomMapModified);
        int i = 0;
        for (Map.Entry<Integer, AtomType> entry : atomIdentifierMapGAFF.entrySet()) {
            String value = entry.getValue().toString();
            // String valueOPLS = atomIdentifierMapOPLS.get(i).toString();

            value = value.replace("AtomType[", "").replace("]", "");
            // valueOPLS = valueOPLS.replace("AtomType[", "").replace("]", "");
            modifiedAtomIdentifierMapGraphene.put(entry.getKey(), value);
            //modifiedAtomIdentifierMapOPLS.put(i, valueOPLS);
        }
        //   System.out.println(modifiedAtomIdentifierMapGraphene + " modified Identifier");
        return modifiedAtomIdentifierMapGraphene;
    }

    Map<String, AtomType> typeMapNew = new HashMap<>();

    protected Map<Integer, AtomType> atomIdentifierGraphene(Map<Integer, String> atomMapModified) {
        Map<Integer, AtomType> atomIdentifierMapGraphene = new HashMap<>();
        int n = 0;
        if (atomMapModified.get(0) == null) {
            n = 1;
        }
        int numHydr = 0, numCarb = 0, numOxy = 0;
        for (int i = n; i < atomMapModified.size(); i++) {
            String atomName = atomMapModified.get(i);
            if (atomName.equals("CX")) {
                AtomType CX = new AtomType(Carbon.INSTANCE, "CX");
                atomIdentifierMapGraphene.put(i, CX);
            } else if (atomName.equals("CY")) {
                AtomType CY = new AtomType(Carbon.INSTANCE, "CY");
                atomIdentifierMapGraphene.put(i, CY);
            } else if (atomName.equals("HK")) {
                AtomType HK = new AtomType(Hydrogen.INSTANCE, "HK");
                atomIdentifierMapGraphene.put(i, HK);
            } else if (atomName.equals("C4")) {
                AtomType C4 = new AtomType(Carbon.INSTANCE, "C4");
                atomIdentifierMapGraphene.put(i, C4);
            } else if (atomName.equals("OJ")) {
                AtomType OJ = new AtomType(Oxygen.INSTANCE, "OJ");
                atomIdentifierMapGraphene.put(i, OJ);
            } else if (atomName.equals("OK")) {
                AtomType OK = new AtomType(Oxygen.INSTANCE, "OK");
                atomIdentifierMapGraphene.put(i, OK);
            } else if (atomName.equals("OL")) {
                AtomType OL = new AtomType(Oxygen.INSTANCE, "OL");
                atomIdentifierMapGraphene.put(i, OL);
            } else if (atomName.equals("OE")) {
                AtomType OE = new AtomType(Oxygen.INSTANCE, "OE");
                atomIdentifierMapGraphene.put(i, OE);
            } else if (atomName.equals("H1")) {
                AtomType H1 = new AtomType(Hydrogen.INSTANCE, "H1");
                atomIdentifierMapGraphene.put(i, H1);
            } else if (atomName.equals("H13")) {
                AtomType H13 = new AtomType(Hydrogen.INSTANCE, "H13");
                atomIdentifierMapGraphene.put(i, H13);
            } else if (atomName.equals("H15")) {
                AtomType H15 = new AtomType(Hydrogen.INSTANCE, "H15");
                atomIdentifierMapGraphene.put(i, H15);
            } else if (atomName.equals("H5")) {
                AtomType H5 = new AtomType(Hydrogen.INSTANCE, "H5");
                atomIdentifierMapGraphene.put(i, H5);
            } else if (atomName.equals("H6")) {
                AtomType H6 = new AtomType(Hydrogen.INSTANCE, "H6");
                atomIdentifierMapGraphene.put(i, H6);
            } else if (atomName.equals("H7")) {
                AtomType H7 = new AtomType(Hydrogen.INSTANCE, "H7");
                atomIdentifierMapGraphene.put(i, H7);
            } else if (atomName.equals("H8")) {
                AtomType H8 = new AtomType(Hydrogen.INSTANCE, "H8");
                atomIdentifierMapGraphene.put(i, H8);
            } else if (atomName.equals("H9")) {
                AtomType H9 = new AtomType(Hydrogen.INSTANCE, "H9");
                atomIdentifierMapGraphene.put(i, H9);
            } else if (atomName.equals("H10")) {
                AtomType H10 = new AtomType(Hydrogen.INSTANCE, "H10");
                atomIdentifierMapGraphene.put(i, H10);
            } else if (atomName.equals("H11")) {
                AtomType H11 = new AtomType(Hydrogen.INSTANCE, "H11");
                atomIdentifierMapGraphene.put(i, H11);
            } else if (atomName.equals("H12")) {
                AtomType H12 = new AtomType(Hydrogen.INSTANCE, "H12");
                atomIdentifierMapGraphene.put(i, H12);
            } else if (atomName.equals("H14")) {
                AtomType H14 = new AtomType(Hydrogen.INSTANCE, "H14");
                atomIdentifierMapGraphene.put(i, H14);
            } else if (atomName.equals("H2")) {
                AtomType H2 = new AtomType(Hydrogen.INSTANCE, "H2");
                atomIdentifierMapGraphene.put(i, H2);
            } else if (atomName.equals("H3")) {
                AtomType H3 = new AtomType(Hydrogen.INSTANCE, "H3");
                atomIdentifierMapGraphene.put(i, H3);
            } else if (atomName.equals("H4")) {
                AtomType H4 = new AtomType(Hydrogen.INSTANCE, "H4");
                atomIdentifierMapGraphene.put(i, H4);
            }

        }
        return atomIdentifierMapGraphene;
    }

    public AtomType returnElement(String elementName) {
        elementReceiverMap.put("CX", new AtomType(Carbon.INSTANCE, "CX"));
        elementReceiverMap.put("CY", new AtomType(Carbon.INSTANCE, "CY"));
        elementReceiverMap.put("CZ", new AtomType(Carbon.INSTANCE, "CZ"));
        elementReceiverMap.put("CA", new AtomType(Carbon.INSTANCE, "CA"));
        elementReceiverMap.put("C4", new AtomType(Carbon.INSTANCE, "C4"));
        elementReceiverMap.put("C", new AtomType(Carbon.INSTANCE, "C"));
        elementReceiverMap.put("HO", new AtomType(Hydrogen.INSTANCE, "HO"));
        elementReceiverMap.put("HK", new AtomType(Hydrogen.INSTANCE, "HK"));
        elementReceiverMap.put("OH", new AtomType(Oxygen.INSTANCE, "OH"));
        elementReceiverMap.put("OJ", new AtomType(Oxygen.INSTANCE, "OJ"));
        elementReceiverMap.put("OK", new AtomType(Oxygen.INSTANCE, "OK"));
        elementReceiverMap.put("O2", new AtomType(Oxygen.INSTANCE, "O2"));
        elementReceiverMap.put("OL", new AtomType(Oxygen.INSTANCE, "OL"));
        elementReceiverMap.put("OS", new AtomType(Oxygen.INSTANCE, "OS"));
        elementReceiverMap.put("O", new AtomType(Oxygen.INSTANCE, "O"));
        elementReceiverMap.put("OE", new AtomType(Oxygen.INSTANCE, "OE"));
        elementReceiverMap.put("N1", new AtomType(Nitrogen.INSTANCE, "N1"));
        elementReceiverMap.put("N2", new AtomType(Nitrogen.INSTANCE, "N2"));
        elementReceiverMap.put("N3", new AtomType(Nitrogen.INSTANCE, "N3"));
        elementReceiverMap.put("H1", new AtomType(Hydrogen.INSTANCE, "H1"));
        elementReceiverMap.put("H2", new AtomType(Hydrogen.INSTANCE, "H2"));
        elementReceiverMap.put("H3", new AtomType(Hydrogen.INSTANCE, "H3"));
        elementReceiverMap.put("H4", new AtomType(Hydrogen.INSTANCE, "H4"));
        elementReceiverMap.put("H5", new AtomType(Hydrogen.INSTANCE, "H5"));
        elementReceiverMap.put("H6", new AtomType(Hydrogen.INSTANCE, "H6"));
        elementReceiverMap.put("H7", new AtomType(Hydrogen.INSTANCE, "H7"));
        elementReceiverMap.put("H8", new AtomType(Hydrogen.INSTANCE, "H8"));
        elementReceiverMap.put("H9", new AtomType(Hydrogen.INSTANCE, "H9"));
        elementReceiverMap.put("H10", new AtomType(Hydrogen.INSTANCE, "H10"));
        elementReceiverMap.put("H11", new AtomType(Hydrogen.INSTANCE, "H11"));
        elementReceiverMap.put("H12", new AtomType(Hydrogen.INSTANCE, "H12"));
        elementReceiverMap.put("H13", new AtomType(Hydrogen.INSTANCE, "H13"));
        elementReceiverMap.put("H14", new AtomType(Hydrogen.INSTANCE, "H14"));
        elementReceiverMap.put("H15", new AtomType(Hydrogen.INSTANCE, "H15"));

        return elementReceiverMap.get(elementName);
    }

    public double[] atomicPotGraphene(String atomType) {
        HashMap<String, double[]> atomicConstant = new HashMap<>();
        atomicConstant.put("CA", new double[]{3.851, 0.105});
        atomicConstant.put("CY", new double[]{3.851, 0.105});
        atomicConstant.put("CX", new double[]{3.851, 0.105});//1.9080 , 0.0860
        atomicConstant.put("C4", new double[]{3.851, 0.105});
        atomicConstant.put("C", new double[]{3.5, 0.105});
        atomicConstant.put("OE", new double[]{3.5, 0.2104});
        atomicConstant.put("O2", new double[]{3.5, 0.2104});
        atomicConstant.put("OJ", new double[]{3.5, 0.2104});
        atomicConstant.put("OK", new double[]{3.5, 0.2104});
        atomicConstant.put("OL", new double[]{3.5, 0.2104});
        atomicConstant.put("OS", new double[]{3.5, 0.1700});
        atomicConstant.put("OH", new double[]{3.5, 0.2104});
        atomicConstant.put("O", new double[]{3.5, 0.2100});
        atomicConstant.put("HK", new double[]{2.886, 0.044});
        atomicConstant.put("H1", new double[]{2.886, 0.044});
        atomicConstant.put("H2", new double[]{2.886, 0.044});
        atomicConstant.put("H3", new double[]{2.886, 0.044});
        atomicConstant.put("H4", new double[]{2.886, 0.044});
        atomicConstant.put("H5", new double[]{2.886, 0.044});
        atomicConstant.put("H6", new double[]{2.886, 0.044});
        atomicConstant.put("H7", new double[]{2.886, 0.044});
        atomicConstant.put("H8", new double[]{2.886, 0.044});
        atomicConstant.put("H9", new double[]{2.886, 0.044});
        atomicConstant.put("H10", new double[]{2.886, 0.044});
        atomicConstant.put("H11", new double[]{2.886, 0.044});
        atomicConstant.put("H12", new double[]{2.886, 0.044});
        atomicConstant.put("H13", new double[]{2.886, 0.044});
        atomicConstant.put("H14", new double[]{2.886, 0.044});
        atomicConstant.put("H15", new double[]{2.886, 0.044});
        atomicConstant.put("N2", new double[]{3.66, 0.069});
        double[] sample = atomicConstant.get(atomType);
        return sample;
    }

    public static void main(String[] args) {
        PDBReaderGraphene pdbReaderGraphene = new PDBReaderGraphene();
        String confName = "F://Avagadro//pristine";
        pdbReaderGraphene.readPDBGrapheneFile(confName);
    }
}
