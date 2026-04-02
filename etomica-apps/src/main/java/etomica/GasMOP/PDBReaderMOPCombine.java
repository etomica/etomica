package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.chem.elements.*;
import etomica.molecule.IMolecule;
import etomica.potential.UFF.PDBReaderMOP;
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

public class PDBReaderMOPCombine {
    Map<Integer,String> atomMap = new HashMap<>();
    Map<Integer,String> atomMapModified = new HashMap<>();
    List<Integer> listFlourideConnect = new ArrayList<>();
    public Map<Integer, String> modifiedAtomIdentifierMap = new HashMap<>();
    Map<Integer,AtomType> atomIdentifierMap = new HashMap<>();
    public Map<Integer,Vector> positions = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectivity = new ArrayList<>();
    public Map<String, AtomType> typeMapNew = new HashMap<>();
    ISpecies species;
    int connect = 0;
    public void readPDBFile(String confName) {
        int m=0;
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
            atomMap.clear();
            Integer atomNumber = 0;
            while ((line = bufReader.readLine()) != null) {
                parseLineReader(line, typeMap, atomMap, positions);
                m++;
                if (line.startsWith("CONECT")) {
                    String[] parts = line.trim().split("\\s+");
                    if (parts[1].length() > 4){
                        if (parts[1].startsWith("9") ){
                            String chunk = parts[1].substring(0,  4);
                            atomNumber = Integer.parseInt(chunk);
                        }else {
                            String chunk = parts[1].substring(0,  5);
                            atomNumber = Integer.parseInt(chunk);
                        }
                    }  else {
                        atomNumber = Integer.parseInt(parts[1]);
                    }

                    if (currentAtomList == null || atomNumber != currentAtomList.get(0)) {
                        // start a new list for the current atom
                        currentAtomList = new ArrayList<>();
                        connectivity.add(currentAtomList);
                        currentAtomList.add(atomNumber);
                    }
                    for (int i = 2; i < parts.length; i++) {
                        String part = parts[i];

                        if (part.startsWith("9")) {
                            if (part.length() <= 4) {
                                currentAtomList.add(Integer.parseInt(part));
                            } else {
                                for (int j = 0; j + 4 <= part.length(); j += 4) {
                                    String chunk = part.substring(j, j + 4);
                                    currentAtomList.add(Integer.parseInt(chunk));
                                }
                            }
                        } else if (part.startsWith("1")) {
                            if (part.length() <= 5) {
                                currentAtomList.add(Integer.parseInt(part));
                            } else {
                                for (int j = 0; j + 5 <= part.length(); j += 5) {
                                    String chunk = part.substring(j, j + 5);
                                    currentAtomList.add(Integer.parseInt(chunk));
                                }
                            }
                        } else {
                            // Fallback for other unexpected cases
                            currentAtomList.add(Integer.parseInt(part));
                        }
                    }
                    connect++;
                }
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        setPositions(positions);
        setConnectivity(connectivity);
        setAtomMap(atomMap);
    }

    public void setAtomMap(Map<Integer, String> atomMap){
        this.atomMap = atomMap;
    }
    public Map<Integer, String> getAtomMap(){return atomMap;}
    public void setAtomMapModified(Map<Integer, String> atomMapModified){
        this.atomMapModified = atomMapModified;
    }
    public Map<Integer, String> getAtomMapModified(){return atomMapModified;}


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
            String symbol = line.substring(12, 14).trim();
            symbol = symbol.replaceAll("\\d{1,2}$", "");
            int atomNumber = Integer.parseInt(line.substring(7,11).trim()) ;
            positions.put(atomNumber, positn);
            //  System.out.println(positions.get(1));
            // System.out.println(atomNumber+" "+positn);
            if (atomMap.containsKey(atomNumber)){
                // System.out.println("Error");
                //  System.out.println("Duplicate atomNumber detected: " + atomNumber);
            } else {
                atomMap.put(atomNumber, symbol);
                // System.out.println("Added atomNumber: " + atomNumber + " with symbol: " + symbol);
            }
        }
    }

    public ISpecies getSpecies (PDBReaderMOP pdbReaderMOP, String confName, boolean isDynamic, Vector centreMOP, boolean setMOPmassInfinite, boolean isMOPVirialAbsent){
        double massSum = 0;
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        SpeciesBuilder speciesBuilderNewMod =  new SpeciesBuilder(Space3D.getInstance());
        ArrayList<Integer> bondList = new ArrayList<>();
        AutoMOPCombines autoMOPCombines = new AutoMOPCombines();
        AtomType typeNew;
        readPDBFile(confName);
        ArrayList<ArrayList<Integer>> connectedAtoms =getConnectivity();
        ArrayList<ArrayList<Integer>> connectivityTemp = getConnectivityTemp(connectedAtoms);
       // PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
       // pdbReaderMOP.setConnectivity(connectedAtoms);
        Map<Integer,String> atomMap = getAtomMap();
        // System.out.println(atomMap + ": atomMap");
        atomMapModified =pdbReaderMOP.getatomMapModified(atomMap);
        Map<Integer, String> atomIdentifierMapMod = atomIdentifierMapModified(connectivityTemp, atomMap);

        ArrayList<ArrayList<Integer>> connectivityModified = pdbReaderMOP.getconnectivityModified(connectivityTemp);
        // System.out.println(connectivityModified+ ": connectedAtomModified" );
        Vector dr = Vector.d(centreMOP.getD());
        for(int i = 0; i < atomIdentifierMapMod.size(); i++) {
            String nameNew = String.valueOf(atomIdentifierMapMod.get(i));
          //  int startIndex = symbol.indexOf("[") + 1;
          //  int endIndex = symbol.indexOf("]");
          //  String nameNew = symbol.substring(startIndex, endIndex);
            // System.out.println(atomIdentifierMap.get(i));
            AtomType newName = pdbReaderMOP.returnElement(nameNew, setMOPmassInfinite);
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
            for(int i=0; i<atomIdentifierMapMod.size(); i++){
                IAtom a = children.get(i);
                String nameNew = String.valueOf(atomIdentifierMapMod.get(i));
             /*   int startIndex = symbol.indexOf("[") + 1;
                int endIndex = symbol.indexOf("]");
                String nameNew = symbol.substring(startIndex, endIndex);*/
                typeNew = typeMapNew.get(nameNew);
                Vector v = a.getPosition();
                v.ME(centreMOP);
                // System.out.println(v);
                speciesBuilderNewMod.addAtom(typeNew,v, "" );
            }
            species = speciesBuilderNewMod.setDynamic(isDynamic).build();
        }
        return species;
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

    public Map<Integer,String> atomIdentifierMapModified(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified) {
        Map<Integer, AtomType> atomIdentifierMap = atomIdentifier(connectivityModified, atomMapModified);
     /*   for(int i=0; i<atomIdentifierMap.size(); i++){
            String atomType = String.valueOf(atomIdentifierMap.get(i));
            String value = atomType.replace("AtomType[", "").replace("]", "");
            if(value.equals("null")){
                System.out.println(i + " " + atomType);
            }
        }*/
        // System.exit(1);
        System.out.println(atomIdentifierMap);
        for (Map.Entry<Integer, AtomType> entry : atomIdentifierMap.entrySet()) {
            String value = entry.getValue().toString();

            value = value.replace("AtomType[", "").replace("]", "");
            modifiedAtomIdentifierMap.put(entry.getKey(), value);
        }
        //  System.out.println(modifiedAtomIdentifierMap + " modified Identifier");
        return modifiedAtomIdentifierMap;
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
        System.out.println(connectivityModified);
        System.out.println(atomMapModified);
      //  System.exit(1);
        //System.out.println(connectivityModified + "connectivityModified");
        System.out.println(connectivityModified.get(connectivityModified.size()-1));
        for (int i = 0; i < connectivityModified.size(); i++) {
            Integer retriveArrayFirstElementNumber = connectivityModified.get(i).get(0);
            //System.out.println(atomMapModified + "atomMap Modified");
            //System.out.println(connectivityModified + " connectivityModified");
            String retriveArrayFirstElementName = atomMapModified.get(retriveArrayFirstElementNumber);
            System.out.println(i +" "+retriveArrayFirstElementName + " :retriveFirstName");
            System.out.println(retriveArrayFirstElementName + "atomMapModified");
            //Carbon atoms
            if (retriveArrayFirstElementName.equals("C")) {
                int arrayListSize = connectivityModified.get(i).size();
                //Alkane checker
                if (arrayListSize == 5) {
                    AtomType C_3 = new AtomType(Carbon.INSTANCE, "C_3");
                    //System.out.println("The atom " + (i) +" is C_3" );
                    atomIdentifierMap.put(i, C_3);
                }

                if (arrayListSize == 4) {
                    AtomType C_2 = new AtomType(Carbon.INSTANCE, "C_2");
                    //System.out.println("The atom " + (i) +" is C_2 " );
                    atomIdentifierMap.put(i, C_2);
                }
                if (arrayListSize == 3 || arrayListSize ==2) {
                    AtomType C_1 = new AtomType(Carbon.INSTANCE, "C_1");
                    //System.out.println("The atom " + (i)+" is C_1 "  );
                    atomIdentifierMap.put(i, C_1);
                }else {
                    AtomType C_3 = new AtomType(Carbon.INSTANCE, "C_3");
                    //System.out.println("The atom " + (i) +" is C_3" );
                    atomIdentifierMap.put(i, C_3);
                }
            } else if (retriveArrayFirstElementName.equals("CL") || retriveArrayFirstElementName.equals("Cl")) {
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
                AtomType H = new AtomType(Hydrogen.INSTANCE, "H");
                //System.out.println("The atom " + (i)+" is H " );
                atomIdentifierMap.put(i, H);
            }
            else if (retriveArrayFirstElementName.equals("O")) {
                int arrayListSize = connectivityModified.get(i).size();
                if(arrayListSize == 2){
                    AtomType O_1 = new AtomType(Oxygen.INSTANCE, "O_1");
                    //System.out.println("The atom " + (i)+" is O_2 " );
                    atomIdentifierMap.put(i, O_1);
                } else if (arrayListSize == 3) {

                    AtomType O_3 =AtomType.element(Oxygen.INSTANCE, "O_3");
                    //System.out.println("The atom " + (i ) + " is O_3 ");
                    atomIdentifierMap.put(i, O_3);
                }  else if (arrayListSize == 4) {

                    AtomType O_2 =AtomType.element(Oxygen.INSTANCE, "O_2");
                    //System.out.println("The atom " + (i ) + " is O_3 ");
                    atomIdentifierMap.put(i, O_2);
                }else {
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
                    AtomType N_1 = new AtomType(Nitrogen.INSTANCE, "N_1");
                    //System.out.println("The atom " + (i) +" is N_1 " );
                    atomIdentifierMap.put(i, N_1);

                } else if (arrayListSize ==4 && oxygenCounter>1) {
                    AtomType N_2 = new AtomType(Nitrogen.INSTANCE, "N_2");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, N_2);
                } else {
                    AtomType N_3 = new AtomType(Nitrogen.INSTANCE, "N_3");
                    //System.out.println("The atom " + (i) +" is N_3 " );
                    atomIdentifierMap.put(i, N_3);
                }

            } else if (retriveArrayFirstElementName.equals("S")) {
                //types are -S-, - -S=O =O, -S=O =O, -S=O =O -OH(-)
                int arrayListSize = connectivityModified.get(i).size();
                if(arrayListSize ==2){
                    AtomType S_2 = new AtomType(Sulfur.INSTANCE, "S_2");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, S_2);
                }else {
                    AtomType S_3 = new AtomType(Sulfur.INSTANCE, "S_3");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, S_3);
                }

            } else if (retriveArrayFirstElementName.equals("P")){
                IElement Phosphorus_3 = null;
                AtomType P_3 = new AtomType(Phosphorus.INSTANCE, "P_3");
                //System.out.println("The atom " + (i+1) +" is P_3 " );
                atomIdentifierMap.put(i, P_3);

            } else if (retriveArrayFirstElementName.equals("SI")) {
                AtomType Si = new AtomType(Silicon.INSTANCE, "Si");
                // System.out.println("The atom " + (i+1) +" is Si " );
                atomIdentifierMap.put(i, Si);

            }  else if (retriveArrayFirstElementName.equals("BR")) {

                AtomType Br = new AtomType(Bromine.INSTANCE, "Br");
                //System.out.println("The atom " + (i) +" is Br " );
                atomIdentifierMap.put(i, Br);

            } else if (retriveArrayFirstElementName.equals("F")) {

                AtomType F = new AtomType(Fluorine.INSTANCE, "F");
                // System.out.println("The atom " + (i+1) +" is F " );
                atomIdentifierMap.put(i, F);

            } else {//Metal ions
                if (retriveArrayFirstElementName.equals("RH")|| retriveArrayFirstElementName.equals("Rh")) {
                    AtomType Rh = new AtomType(Rhodium.INSTANCE, "Rh");
                    // System.out.println("The atom " + (i+1) +" is Rh " );
                    atomIdentifierMap.put(i, Rh);

                }  if (retriveArrayFirstElementName.equals("B")) {
                    AtomType B = new AtomType(Boron.INSTANCE, "B");
                    // System.out.println("The atom " + (i+1) +" is Rh " );
                    atomIdentifierMap.put(i, B);

                } else if (retriveArrayFirstElementName.equals("RU")|| retriveArrayFirstElementName.equals("Ru")) {
                    AtomType Ru = new AtomType(Ruthenium.INSTANCE, "Ru");
                    // System.out.println("The atom " + (i+1) +" is Ru " );
                    atomIdentifierMap.put(i, Ru);

                } else if (retriveArrayFirstElementName.equals("NI") || retriveArrayFirstElementName.equals("Ni")) {
                    AtomType Ni = new AtomType(Nickel.INSTANCE, "Ni");
                    // System.out.println("The atom " + (i+1) +" is Ni " );
                    atomIdentifierMap.put(i, Ni);

                } else if (retriveArrayFirstElementName.equals("CU")|| retriveArrayFirstElementName.equals("Cu")) {
                    AtomType Cu = new AtomType(Copper.INSTANCE, "Cu");
                    //System.out.println("The atom " + (i+1) +" is Cu " );
                    atomIdentifierMap.put(i, Cu);

                } else if (retriveArrayFirstElementName.equals("Fe") || retriveArrayFirstElementName.equals("FE")) { // Not this
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

                } else if (retriveArrayFirstElementName.equals("CO") || retriveArrayFirstElementName.equals("Co")) {
                    AtomType Co = new AtomType(Cobalt.INSTANCE, "Co");
                    //System.out.println("The atom " + (i+1) +" is Co " );
                    atomIdentifierMap.put(i, Co);

                } else if (retriveArrayFirstElementName.equals("CR")|| retriveArrayFirstElementName.equals("Cr")) {
                    AtomType Cr = new AtomType(Chromium.INSTANCE, "Cr");
                    //System.out.println("The atom " + (i+1) +" is Cr " );
                    atomIdentifierMap.put(i, Cr);

                } else if (retriveArrayFirstElementName.equals("PD")|| retriveArrayFirstElementName.equals("Pd")) {
                    AtomType Pd = new AtomType(Palladium.INSTANCE, "Pd");
                    //System.out.println("The atom " + (i+1) +" is Pd " );
                    atomIdentifierMap.put(i, Pd);

                } else if (retriveArrayFirstElementName.equals("MO")|| retriveArrayFirstElementName.equals("Mo")) {
                    AtomType Mo = new AtomType(Molybdenum.INSTANCE, "Mo");
                    // System.out.println("The atom " + (i+1) +" is Mo " );
                    atomIdentifierMap.put(i, Mo);

                }else if (retriveArrayFirstElementName.equals("MN")|| retriveArrayFirstElementName.equals("Mn")) {
                    AtomType Mn = new AtomType(Manganese.INSTANCE, "Mn");
                    // System.out.println("The atom " + (i+1) +" is Mo " );
                    atomIdentifierMap.put(i, Mn);

                } else if (retriveArrayFirstElementName.equals("ZR")|| retriveArrayFirstElementName.equals("Zr")) {
                    AtomType Zr = new AtomType(Zirconium.INSTANCE, "Zr");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, Zr);

                } else if (retriveArrayFirstElementName.equals("V")) {
                    AtomType V = new AtomType(Vanadium.INSTANCE, "V");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, V);

                } else if (retriveArrayFirstElementName.equals("W")) {// Not this
                    AtomType W = new AtomType(Tungsten.INSTANCE, "W");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, W);

                } else if (retriveArrayFirstElementName.equals("MG")|| retriveArrayFirstElementName.equals("Mg")) {
                    AtomType Mg = new AtomType(Magnesium.INSTANCE, "Mg");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, Mg);

                } else if (retriveArrayFirstElementName.equals("ZN")|| retriveArrayFirstElementName.equals("Zn")) {
                    AtomType Zn = new AtomType(Zinc.INSTANCE, "Zn");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, Zn);
                } else if (retriveArrayFirstElementName.equals("IR")|| retriveArrayFirstElementName.equals("Ir")) {
                    AtomType Ir = new AtomType(Vanadium.INSTANCE, "Ir");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, Ir);
                }  else if (retriveArrayFirstElementName.equals("TI")|| retriveArrayFirstElementName.equals("Ti")) {
                    AtomType Ti = new AtomType(Titanium.INSTANCE, "Ti");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, Ti);
                }  else if (retriveArrayFirstElementName.equals("In")|| retriveArrayFirstElementName.equals("IN")) {
                    AtomType In = new AtomType(Indium.INSTANCE, "In");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, In);
                }  else if (retriveArrayFirstElementName.equals("Ca")|| retriveArrayFirstElementName.equals("CA")) {
                    AtomType Ca = new AtomType(Carbon.INSTANCE, "Ca");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, Ca);
                }  else if (retriveArrayFirstElementName.equals("CD")|| retriveArrayFirstElementName.equals("Cd")) {
                    AtomType Cd = new AtomType(Cadmium.INSTANCE, "Cd");
                    //System.out.println("The atom " + (i) +" is N_2 " );
                    atomIdentifierMap.put(i, Cd);
                } else {
                    System.out.println(atomIdentifierMap.get(retriveArrayFirstElementNumber));
                    throw new RuntimeException("Atom not found " + retriveArrayFirstElementName + " " +retriveArrayFirstElementNumber);
                }
            }
        }
        // System.out.println(atomIdentifierMap + "atomIdentifier");
        return atomIdentifierMap;
    }
    public List<Integer> settingUpFluorideList(ArrayList<ArrayList<Integer>> connectivityMap, Map<Integer, String> atomMap, ISpecies speciesOne){
        List<Integer> fluorideList = new ArrayList<>();
        for (int i =0 ; i < atomMap.size(); i++){
            String atomName = atomMap.get(i);
            if (atomName.equals("F")){
                fluorideList.add(i);
            }
        }
        setListFlourideConnect(fluorideList);
        return fluorideList;
    }
    private void setPositions(Map<Integer, Vector> positions) {
        this.positions = positions;
    }
    public Map<Integer, Vector> getPositions (){
        return positions;
    }
    public void setConnectivity(ArrayList<ArrayList<Integer>>connectivity){this.connectivity = connectivity;}
    public ArrayList<ArrayList<Integer>> getConnectivity(){
        return connectivity;
    }
    public void setListFlourideConnect (List<Integer> connectFluoride){ this.listFlourideConnect = connectFluoride;}
    public List<Integer> getSpeciesConnectList(){return listFlourideConnect;}
}
