package etomica.potential.GAFF;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.chem.elements.*;
import etomica.molecule.IMolecule;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class COFReader {
    public Map<String ,Vector> positions = new HashMap<>();
    public Map<Integer, String> atomMap = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectivity = new ArrayList<>();
    public Map<Integer, String> modifiedAtomIdentifierMap = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectivityModified = new ArrayList<>();
    public Map<Integer, String> atomMapModified = new HashMap<>();
    public Map<Integer, AtomType> atomIdentifierMap = new HashMap<>();
    public Map<String, AtomType> typeMapNew = new HashMap<>();
    public Map<String, AtomType> elementReceiverMap = new HashMap<>();
    int connect = 0;
    public ISpecies species;
    public void readCOFFile(String confName){
        String fileName = confName+".cif";
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

                if (line.startsWith("C") || line.startsWith("N") || line.startsWith("O") || line.startsWith("H")) {
                    String[] parts = line.trim().split("\\s+");
                    int atomNumber = Integer.parseInt(parts[1]);
                    if (currentAtomList == null || atomNumber != currentAtomList.get(0)) {
                        // start a new list for the current atom
                        currentAtomList = new ArrayList<>();
                        connectivity.add(currentAtomList);
                        currentAtomList.add(atomNumber);
                    }
                    for (int i = 2; i < parts.length; i++) {
                        currentAtomList.add(Integer.parseInt(parts[i]));
                    }
                    connect++;
                }
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        setPositions(positions);
    }

    protected void parseLineReader(String line, Map<String, AtomType> typeMap, Map<Integer, String> atomMap, Map<String, Vector> positions){
        line = line.trim();
        if (line.length() < 6) {
            return;
        }
        if (line.startsWith("C") || line.startsWith("N") || line.startsWith("O") || line.startsWith("H")) {
            //coordinates of atom and create atom
            String[] parts = line.split(" ");
            double x = Double.parseDouble(line.substring(13, 20));
            double y = Double.parseDouble(line.substring(24, 29));
            double z = Double.parseDouble(line.substring(33, 39));
            String partO = parts[0];
            Vector positn = Vector.of(x, y, z);
            String symbol = line.substring(12, 16).trim();
            positions.put(partO, positn);
           /* if (atomMap.containsKey(atomNumber)){
                // System.out.println("Error");
            } else {
                atomMap.put(atomNumber, symbol);
            }*/
        }
    }
    private void setPositions(Map<String, Vector> positions) {
        this.positions = positions;
    }

    public ISpecies getSpeciesCOF(String confName, boolean isDynamic, Vector centreCOF, boolean setMassInfinite){
        double massSum = 0;
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        SpeciesBuilder speciesBuilderNewMod =  new SpeciesBuilder(Space3D.getInstance());
        AtomType typeNew;
        readCOFFile(confName);
        ArrayList<ArrayList<Integer>> connectedAtoms =getConnectivity();
        ArrayList<ArrayList<Integer>> connectivityTemp = getConnectivityTemp(connectedAtoms);
        Map<Integer, String> atomIdentifierMapMod = atomIdentifierMapModified(connectivityTemp, atomMap);
        ArrayList<ArrayList<Integer>> connectivityModified = getconnectivityModified(connectivityTemp);
        Map<Integer,String> atomMap = getAtomMap(connectivityTemp);
        atomMapModified =getatomMapModified(atomMap);
        Vector dr = Vector.d(centreCOF.getD());
        for(int i = 0; i < atomIdentifierMap.size(); i++) {
            String symbol = String.valueOf(atomIdentifierMap.get(i));
            int startIndex = symbol.indexOf("[") + 1;
            int endIndex = symbol.indexOf("]");
            String nameNew = symbol.substring(startIndex, endIndex);
            AtomType newName = returnElement(nameNew, setMassInfinite);
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
        if(!setMassInfinite){
            IMolecule molecule = species.makeMolecule();
            IAtomList children = molecule.getChildList();
            int nAtoms = children.size();
            for (int i = 0; i < nAtoms; i++) {
                IAtom a = children.get(i);
                double mass = a.getType().getMass();
                if (massSum == 0) {
                    centreCOF.PEa1Tv1(mass, a.getPosition());
                } else {
                    dr.E(a.getPosition());
                    centreCOF.PEa1Tv1(mass, dr);
                }
                massSum += mass;
            }
            centreCOF.TE(1.0 / massSum);
            for(int i=0; i<atomIdentifierMap.size(); i++){
                IAtom a = children.get(i);
                String symbol = String.valueOf(atomIdentifierMap.get(i));
                int startIndex = symbol.indexOf("[") + 1;
                int endIndex = symbol.indexOf("]");
                String nameNew = symbol.substring(startIndex, endIndex);
                typeNew = typeMapNew.get(nameNew);
                Vector v = a.getPosition();
                v.ME(centreCOF);
                speciesBuilderNewMod.addAtom(typeNew,v, "" );
            }
            species = speciesBuilderNewMod.setDynamic(isDynamic).build();
        }
        return species;
    }

    public ArrayList<ArrayList<Integer>> getConnectivity(){
        return connectivity;
    }

    public ArrayList<ArrayList<Integer>> getConnectivityTemp(ArrayList<ArrayList<Integer>>connectivity){
        ArrayList<ArrayList<Integer>> connectivityTemp = new ArrayList<>();
        for(int i =0; i<connectivity.size(); i++){
            ArrayList<Integer> uniqueList = new ArrayList<>();
            for (Integer number : connectivity.get(i)) {
                if (!uniqueList.contains(number)) {
                    uniqueList.add(number);
                }
            }
            connectivityTemp.add(i, uniqueList);
        }
        return connectivityTemp;
    }

    public Map<Integer,String> atomIdentifierMapModified (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified) {
        Map<Integer, AtomType> atomIdentifierMap = atomIdentifier(connectivityModified, atomMapModified);

        for (Map.Entry<Integer, AtomType> entry : atomIdentifierMap.entrySet()) {
            String value = entry.getValue().toString();

            value = value.replace("AtomType[", "").replace("]", "");
            modifiedAtomIdentifierMap.put(entry.getKey(), value);
        }
        //  System.out.println(modifiedAtomIdentifierMap + " modified Identifier");
        return modifiedAtomIdentifierMap;
    }

    protected Map<Integer, AtomType> atomIdentifier(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        Map<Integer, AtomType> atomIdentifierMap = new HashMap<>();

        return atomIdentifierMap;
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

    public Map<Integer, String> getatomMapModified( Map<Integer, String> atomMap){
        for (Map.Entry<Integer, String> entry : atomMap.entrySet()) {
            atomMapModified.put(entry.getKey() - 1, entry.getValue());
        }
        return atomMapModified;
    }

    public AtomType returnElement(String elementName, boolean isInfinite){
        if(isInfinite){
            switch (elementName){
                case "H":
                    AtomType H = new AtomType(new ElementSimple("H", Double.POSITIVE_INFINITY), "H");
                    elementReceiverMap.put("H", H);
                    break;
                case "Cu":
                    AtomType Cu = new AtomType(new ElementSimple("Cu", Double.POSITIVE_INFINITY), "Cu");
                    elementReceiverMap.put("Cu", Cu);
                    break;
                case "C_3":
                    AtomType C_3 = new AtomType(new ElementSimple("C_3", Double.POSITIVE_INFINITY), "C_3");
                    elementReceiverMap.put("C_3", C_3);
                    break;
                case "C_2":
                    AtomType C_2 = new AtomType(new ElementSimple("C_2", Double.POSITIVE_INFINITY), "C_2");
                    elementReceiverMap.put("C_2", C_2);
                    break;
                case "C_1":
                    AtomType C_1 = new AtomType(new ElementSimple("C_1", Double.POSITIVE_INFINITY), "C_1");
                    elementReceiverMap.put("C_1", C_1);
                    break;
                case "O_3":
                    AtomType O_3 = new AtomType(new ElementSimple("O_3", Double.POSITIVE_INFINITY), "O_3");
                    elementReceiverMap.put("O_3", O_3);
                    break;
                case "O_2":
                    AtomType O_2 = new AtomType(new ElementSimple("O_2", Double.POSITIVE_INFINITY), "O_2");
                    elementReceiverMap.put("O_2", O_2);
                    break;
                case "O_1":
                    AtomType O_1 = new AtomType(new ElementSimple("O_1", Double.POSITIVE_INFINITY), "O_1");
                    elementReceiverMap.put("O_1", O_1);
                    break;
                case "Ar":
                    AtomType Ar = new AtomType(new ElementSimple("Ar", Double.POSITIVE_INFINITY), "Ar");
                    elementReceiverMap.put("Ar", Ar);
                    break;
            }
          /*/  elementReceiverMap.put("H", AtomType.simple(String.valueOf(new ElementSimple("H")), Double.POSITIVE_INFINITY));
            elementReceiverMap.put("Cu", AtomType.simple(String.valueOf(new ElementSimple("Cu")), Double.POSITIVE_INFINITY));
            elementReceiverMap.put("C_3", AtomType.simple(String.valueOf(new ElementSimple("C_3")), Double.POSITIVE_INFINITY));
            elementReceiverMap.put("C_2", AtomType.simple(String.valueOf(new ElementSimple("C_2")), Double.POSITIVE_INFINITY));
            elementReceiverMap.put("C_1", AtomType.simple(String.valueOf(new ElementSimple("C_1")), Double.POSITIVE_INFINITY));
            elementReceiverMap.put("O_3", AtomType.simple(String.valueOf(new ElementSimple("O_3")), Double.POSITIVE_INFINITY));
            elementReceiverMap.put("O_2", AtomType.simple(String.valueOf(new ElementSimple("O_2")), Double.POSITIVE_INFINITY));
            elementReceiverMap.put("O_1", AtomType.simple(String.valueOf(new ElementSimple("O_1")), Double.POSITIVE_INFINITY));*/
        }else {
            elementReceiverMap.put("C_3", new AtomType(Carbon.INSTANCE, "C_3"));
            elementReceiverMap.put("H", new AtomType(Carbon.INSTANCE, "H"));
            elementReceiverMap.put("C_2", new AtomType(Carbon.INSTANCE, "C_2"));
            elementReceiverMap.put("C_Ar", new AtomType(Carbon.INSTANCE, "C_Ar"));
            elementReceiverMap.put("C_1", new AtomType(Carbon.INSTANCE, "C_1"));
            elementReceiverMap.put("O_3", new AtomType(Oxygen.INSTANCE, "O_3"));
            elementReceiverMap.put("O_1", new AtomType(Oxygen.INSTANCE, "O_1"));
            elementReceiverMap.put("O_Ar", new AtomType(Oxygen.INSTANCE, "O_Ar"));
            elementReceiverMap.put("O_2", new AtomType(Oxygen.INSTANCE, "O_2"));
            elementReceiverMap.put("N_Ar", new AtomType(Nitrogen.INSTANCE, "N_Ar"));
            elementReceiverMap.put("N_3", new AtomType(Nitrogen.INSTANCE, "N_3"));
            elementReceiverMap.put("N_2", new AtomType(Nitrogen.INSTANCE, "N_2"));
            elementReceiverMap.put("N_1", new AtomType(Nitrogen.INSTANCE, "N_1"));
            elementReceiverMap.put("S_3", new AtomType(Sulfur.INSTANCE, "S_3"));
            elementReceiverMap.put("S_2", new AtomType(Sulfur.INSTANCE, "S_2"));
            elementReceiverMap.put("Cl", new AtomType(Chlorine.INSTANCE, "Cl"));
            elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));
            elementReceiverMap.put("Co", new AtomType(Cobalt.INSTANCE, "Co"));
            elementReceiverMap.put("Fe", new AtomType(Iron.INSTANCE, "Fe"));
            elementReceiverMap.put("Ni", new AtomType(Nickel.INSTANCE, "Ni"));
            elementReceiverMap.put("Rh", new AtomType(Rhodium.INSTANCE, "Rh"));
            elementReceiverMap.put("Ru", new AtomType(Ruthenium.INSTANCE, "Ru"));
            elementReceiverMap.put("Ti", new AtomType(Titanium.INSTANCE, "Ti"));
            elementReceiverMap.put("Sc", new AtomType(Scandium.INSTANCE, "Sc"));
            elementReceiverMap.put("Zn", new AtomType(Zinc.INSTANCE, "Zn"));
           /* elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));
            elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));
            elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));
            elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));
            elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));
            elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));
            elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));
            elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));
            elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));
            elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));
            elementReceiverMap.put("V", new AtomType(Vanadium.INSTANCE, "V"));*/
            elementReceiverMap.put("Zr", new AtomType(Zirconium.INSTANCE, "Zr"));
            elementReceiverMap.put("Br", new AtomType(Bromine.INSTANCE, "Br"));
            elementReceiverMap.put("Cr", new AtomType(Chromium.INSTANCE, "Cr"));
            elementReceiverMap.put("Mn", new AtomType(Manganese.INSTANCE, "Mn"));
            elementReceiverMap.put("Cu", new AtomType(Copper.INSTANCE, "Cu"));
            elementReceiverMap.put("Pd", new AtomType(Palladium.INSTANCE, "Pd"));
            elementReceiverMap.put("Ar", new AtomType(Argon.INSTANCE,"Ar"));
            elementReceiverMap.put("He", new AtomType(Helium.INSTANCE,"He"));
            elementReceiverMap.put("Fe_3", new AtomType(Iron.INSTANCE,"Fe_3"));
        }

        return elementReceiverMap.get(elementName);
    }


    public static void main(String[] args) {
        COFReader cofReader = new COFReader();
        Vector postn = new Vector3D(0,0,0);
        String confName = "D://SemVII//COFs//COF/v1/2D-NiPc-BTDA-COF";
        ISpecies species = cofReader.getSpeciesCOF(confName, true, postn, false );
    }
}
