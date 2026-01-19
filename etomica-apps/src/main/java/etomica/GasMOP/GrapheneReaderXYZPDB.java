package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.potential.UFF.GeneralGrapheneReader;
import etomica.potential.UFF.PDBReaderMOP;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.util.collections.IntArrayList;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class GrapheneReaderXYZPDB {
    private ISpecies species;
    double massSum;
    public Map<Integer, String> atomMapGraphene = new HashMap<>();
    public Map<Integer,Vector> positions = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectivityGrapehene = new ArrayList<>();
    public Map<Integer, String> modifiedAtomIdentifierMapGraphene = new HashMap<>();
    public Map<Integer, Vector> positionsGraphene = new HashMap<>();
    private Map<Integer,List<int[]>> bondTypeMap = new HashMap<>();
    private Map<Integer,List<int[]>> angleTypeMap = new HashMap<>();
    private Map<Integer,List<int[]>> torsionTypeMap = new HashMap<>();
    private Map<Integer,List<int[]>> diheadralTypeMap = new HashMap<>();
    private Map<Integer,List<int[]>> improperTypeMap = new HashMap<>();
    public Map<String, AtomType> typeMap = new HashMap<>();
    public List<int[]> dupletsSorted, tripletsSorted, quadrupletsSorted = new ArrayList<>();
    Map<String, AtomType> typeMapNew = new HashMap<>();
    protected Map<Integer, String> atomMap = new HashMap<>();
    protected Set<String> uniqueElements = new HashSet<>();
    protected Map<Integer, String> uniqueAtomTypeMap;
    protected Map<Integer,double[]> bondPotential = new HashMap<>();
    protected Map<Integer,double[]> anglePotential = new HashMap<>();
    protected Map<Integer,double[]> torsionPotential = new HashMap<>();
    protected Map<Integer, Integer> elementNumMap = new HashMap<>();
    protected Map<Integer, Double> chargeMap = new HashMap<>();
    protected Map<Integer, Vector> positionMap = new HashMap<>();
    double Q_TOL = 0.0001;
    protected IntArrayList[] hydroxyIntArr, epoxyIntArr, carboxyIntArr;

    Map<Integer, String> atomTypeMap;
    Map<Integer, String> expandedAtomTypeMap = new HashMap<>();
    Map<String, Double> advChargeMap = new HashMap<>();

    /*public ISpecies getSpecies(String confName, Vector vectorCOM, boolean isInfinite){
        AtomType typeNew;
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        PDBReaderMOP pdbReader = new PDBReaderMOP();
        GrapheneReaderXYZPDB grapheneReader = new GrapheneReaderXYZPDB();
        grapheneReader.readDatafile(confName);
        positionsGraphene = pdbReader.getPositions();
        atomMap = pdbReader.getAtomMap();
        Space space = Space3D.getInstance();
        Vector center = space.makeVector();
        Vector dr = Vector.d(center.getD());
        int n =0;
        for(int i = n; i < atomMap.size(); i++) {
            String symbol = String.valueOf(atomMap.get(i));
            //System.out.println(symbol);
            // System.out.println(atomIdentifierMap.get(i));
            AtomType newName = generalGrapheneReader.returnElement(symbol, isInfinite);
            if (typeMapNew.containsKey(newName.getName())) {
                typeNew = typeMapNew.get(newName.getName());
            } else {
                typeNew = newName;
                typeMapNew.put(newName.getName(), typeNew);
            }

            Vector position = positionsGraphene.get(i);
            //System.out.println(position +  " " + typeNew + "  " + i);
            speciesBuilderNew.addAtom(typeNew, position,  "");
        }
        //System.out.println(typeMapNew + " typeMapNew");
        species= speciesBuilderNew.setDynamic(true).build();
        //   System.out.println(species.getMass() + " first");
        return species;
    }*/

    public void bondingGraphene(String fileName){
        ArrayList<ArrayList<Integer>> connectivity = new ArrayList<>();
        Map<Integer, List<int[]>> bondTypeMap = new HashMap<>();
        Map<Integer, List<int[]>> angleTypeMap = new HashMap<>();
        Map<Integer, List<int[]>> torsionTypeMap = new HashMap<>();
    }




    // Put this enum somewhere in your class (or as a private static enum)

    // Replace your readDatafile(...) with this.
// Notes:
// 1) It reads confName + ".data" (change extension if needed).
// 2) It stores bond/angle coeffs into maps as soon as it is inside that section.
// 3) It still calls your existing parseLineReader(...) for atom/xyz-style lines if you want.
// 4) Add/modify the calls in the switch blocks to match your own data structures.
    public void readDatafile(String confName) {

        // ---- choose correct extension ----
        String fileName = confName + ".data";   // <-- you said .data file (not .xyz)
        // If you really want: confName+".xyz" then change back.

        // ---- your maps (use your existing ones if already fields) ----
        // If these are fields already, DELETE these local declarations.
        Map<Integer, double[]> bondPotentials  = new HashMap<>();
        Map<Integer, double[]> anglePotentials = new HashMap<>();
        Map<Integer, double[]> torsionPotentials = new HashMap<>();
        Map<Integer, double[]> improperPotentials = new HashMap<>();
        Map<Integer, double[]> dihedralPotentials = new HashMap<>();
        List<int[]> bonds = new ArrayList<>();
        int[] bondsArray = new int[2];
        int[] anglesArray = new int[3];
        int[] torsionsArray = new int[4];
        Section section = Section.NONE;
        boolean inMasses = false;
        Map<Integer, Double> massByTypeId = new LinkedHashMap<>(); // keeps file order
        try (BufferedReader bufReader = new BufferedReader(new FileReader(fileName))) {
            String raw;

            while ((raw = bufReader.readLine()) != null) {
                String line = raw.trim();

                // -------------------------
                // skip blanks / comments
                // -------------------------
                if (line.isEmpty() || line.startsWith("#")) continue;


                if (line.equalsIgnoreCase("Masses")) {
                    inMasses = true;
                    continue;
                }

                // if we were in Masses and we hit another section header, stop reading Masses
                // (LAMMPS sections are single words like "Atoms", "Bonds", "Angles", etc.)
                if (inMasses) {
                    String[] t = line.split("\\s+");

                    // If it's not at least 2 tokens, or first token not int, Masses probably ended
                    if (t.length < 2 || !t[0].matches("\\d+")) {
                        inMasses = false;   // fall through to other parsing (or continue)
                        // OPTIONAL: continue; // if you want to ignore this line here
                    } else {
                        int typeId = Integer.parseInt(t[0]);
                        double mass = Double.parseDouble(t[1]);
                        massByTypeId.put(typeId, mass);
                        continue; // done with this line
                    }
                }

                // -------------------------
                // detect section headers
                // -------------------------
                // Coeff sections
                if (equalsHeader(line, "Bond Coeffs"))     { section = Section.bond_Coeffs;    continue; }
                if (equalsHeader(line, "Angle Coeffs"))    { section = Section.angle_Coeffs;   continue; }
                if (equalsHeader(line, "Torsion Coeffs"))  { section = Section.torsion_Coeffs; continue; }
                if (equalsHeader(line, "Dihedral Coeffs")) { section = Section.dihedral_Coeffs;continue; }
                if (equalsHeader(line, "Improper Coeffs")) { section = Section.improper_Coeffs;continue; }

                // Topology sections (optional if your .data has them)
                if (equalsHeader(line, "Atoms"))     { section = Section.atoms;     continue; }
                if (equalsHeader(line, "Bonds"))     { section = Section.bonds;     continue; }
                if (equalsHeader(line, "Angles"))    { section = Section.angles;    continue; }
               // if (equalsHeader(line, "Torsions"))  { section = Section.torsions;  continue; }
                if (equalsHeader(line, "Dihedrals")) { section = Section.dihedrals; continue; }
                if (equalsHeader(line, "Impropers")) { section = Section.impropers; continue; }

                // If line is another random heading (letters but not numeric), exit section safely
                if (looksLikeHeading(line)) {
                    section = Section.NONE;
                    continue;
                }

                // -------------------------
                // parse based on section
                // -------------------------

                switch (section) {

                    case bond_Coeffs: {
                        String[] t = line.split("\\s+");
                        if (t.length >= 3) {
                            int id = Integer.parseInt(t[0]);
                            double k  = Double.parseDouble(t[1]);
                            double r0 = Double.parseDouble(t[2]);
                            bondPotentials.put(id, new double[]{k, r0});
                        }
                        break;
                    }

                    case angle_Coeffs: {
                        String[] t = line.split("\\s+");
                        if (t.length >= 3) {
                            int id = Integer.parseInt(t[0]);
                            double k = Double.parseDouble(t[1]);
                            double theta0 = Double.parseDouble(t[2]);
                            anglePotentials.put(id, new double[]{k, theta0});
                        }
                        break;
                    }

                    case dihedral_Coeffs: {
                        String[] t = line.split("\\s+");
                        if (t.length >= 2) {
                            int id = Integer.parseInt(t[0]);
                            double[] params = new double[t.length - 1];
                            for (int i = 1; i < t.length; i++) {
                                params[i - 1] = Double.parseDouble(t[i]);
                            }
                            torsionPotentials.put(id, params);
                        }
                        break;
                    }

                    case improper_Coeffs: {
                        String[] t = line.split("\\s+");
                        if (t.length >= 1) {
                            int id = Integer.parseInt(t[0]);
                            double[] params = new double[t.length - 1];
                            for (int i = 1; i < t.length; i++) {
                                params[i - 1] = Double.parseDouble(t[i]);
                            }
                            torsionPotentials.put(id, params);
                        }
                        break;
                    }

                    case atoms: {
                        String[] t = line.split("\\s+");

                        // atomID   moleculeID  atomTypeID  partialCharge   X   Y   Z
                        if (t.length >= 5){
                            int numZero = Integer.parseInt(t[0]);
                            int numTwo = Integer.parseInt(t[2]);
                            double numThree = Double.parseDouble(t[3]);
                            double numFour = Double.parseDouble(t[4]);
                            double numFive = Double.parseDouble(t[5]);
                            double numSix = Double.parseDouble(t[6]);
                            positionMap.put(numZero, Vector.of(numFour, numFive, numSix));
                            elementNumMap.put(numZero, numTwo);
                            chargeMap.put(numZero, numThree);
                        }else {
                            throw new RuntimeException("t "+ Arrays.toString(t));
                        }
                        break;
                    }

                    case bonds: {
                        String[] t = line.split("\\s+");
                        if (t.length >= 3) {
                            int bondType = Integer.parseInt(t[1]);
                            int atom1    = Integer.parseInt(t[2]);
                            int atom2    = Integer.parseInt(t[3]);

                            int[] bondPair = new int[]{atom1, atom2};

                            // create list if this bond type is seen first time
                            List<int[]> bondList = bondTypeMap.get(bondType);
                            if (bondList == null) {
                                bondList = new ArrayList<>();
                                bondList.add(bondPair);
                                bondTypeMap.put(bondType, bondList);
                            }else {
                                bondList.add(bondPair);
                            }
                        }
                        break;
                    }

                    case angles: {
                        String[] t = line.split("\\s+");
                        if (t.length >= 4) {
                            int bondType = Integer.parseInt(t[1]);
                            int atom1    = Integer.parseInt(t[2]);
                            int atom2    = Integer.parseInt(t[3]);
                            int atom3    = Integer.parseInt(t[4]);

                            int[] bondPair = new int[]{atom1, atom2, atom3};

                            // create list if this bond type is seen first time
                            List<int[]> bondList = angleTypeMap.get(bondType);
                            if (bondList == null) {
                                bondList = new ArrayList<>();
                                bondList.add(bondPair);
                                angleTypeMap.put(bondType, bondList);
                            }else {
                                bondList.add(bondPair);
                            }
                        }
                        break;
                    }

                    case torsions:{
                        String[] t = line.split("\\s+");
                        if (t.length >= 4) {
                            int bondType = Integer.parseInt(t[1]);
                            int atom1    = Integer.parseInt(t[2]);
                            int atom2    = Integer.parseInt(t[3]);
                            int atom3    = Integer.parseInt(t[4]);
                            int atom4    = Integer.parseInt(t[5]);
                            int[] bondPair = new int[]{atom1, atom2, atom3, atom4};

                            // create list if this bond type is seen first time
                            List<int[]> bondList = torsionTypeMap.get(bondType);
                            if (bondList == null) {
                                bondList = new ArrayList<>();
                                bondList.add(bondPair);
                                torsionTypeMap.put(bondType, bondList);
                            }else {
                                bondList.add(bondPair);
                            }
                        }
                        break;
                    }

                    case dihedrals:{
                        String[] t = line.split("\\s+");
                        if (t.length >= 4) {
                            int bondType = Integer.parseInt(t[1]);
                            int atom1    = Integer.parseInt(t[2]);
                            int atom2    = Integer.parseInt(t[3]);
                            int atom3    = Integer.parseInt(t[4]);
                            int atom4    = Integer.parseInt(t[5]);
                            int[] bondPair = new int[]{atom1, atom2, atom3, atom4};

                            // create list if this bond type is seen first time
                            List<int[]> bondList = diheadralTypeMap.get(bondType);
                            if (bondList == null) {
                                bondList = new ArrayList<>();
                                bondList.add(bondPair);
                                diheadralTypeMap.put(bondType, bondList);
                            }else {
                                bondList.add(bondPair);
                            }
                        }
                        break;
                    }

                    case impropers: { String[] t = line.split("\\s+");
                        if (t.length >= 4) {
                            int bondType = Integer.parseInt(t[1]);
                            int atom1    = Integer.parseInt(t[2]);
                            int atom2    = Integer.parseInt(t[3]);
                            int atom3    = Integer.parseInt(t[4]);
                            int atom4    = Integer.parseInt(t[5]);
                            int[] bondPair = new int[]{atom1, atom2, atom3, atom4};

                            // create list if this bond type is seen first time
                            List<int[]> bondList = improperTypeMap.get(bondType);
                            if (bondList == null) {
                                bondList = new ArrayList<>();
                                bondList.add(bondPair);
                                improperTypeMap.put(bondType, bondList);
                            }else {
                                bondList.add(bondPair);
                            }
                        }
                        break;
                    }
                    default:
                        break;
                }
            }

            setBondTypeMap(bondTypeMap);
            setAngleTypeMap(angleTypeMap);
            setTorsionTypeMap(torsionTypeMap);
            setDihedralTypeMap(diheadralTypeMap);
            setImproperTypeMap(improperTypeMap);
            setConnectivity(connectivityGrapehene);
            Map<Integer, String> uniqueTypes = makeUniqueAtomTypes(massByTypeId);
            setUniqueAtomTypeMap(uniqueTypes);
            atomMap = makeAtomMap(uniqueTypes, elementNumMap);
            setPotentials(bondPotentials, anglePotentials, torsionPotentials, improperPotentials, dihedralPotentials, elementNumMap, chargeMap, atomMap, positionMap);
            System.out.println("Done Reading: " + fileName);

        } catch (IOException e) {
            throw new RuntimeException("Problem reading from " + fileName + ", caught IOException: " + e.getMessage(), e);
        }
    }

    public void chargeGraphene(Map<Integer, String> atomTypeMap, Map<Integer, Double> chargeMap){
        Map<String, List<Integer>> atomsByBaseType = new HashMap<>();

        for (int atomId : atomTypeMap.keySet()) {
            String baseType = atomTypeMap.get(atomId); // e.g. "N_2"
            atomsByBaseType
                    .computeIfAbsent(baseType, k -> new ArrayList<>())
                    .add(atomId);
        }


        for (Map.Entry<String, List<Integer>> entry : atomsByBaseType.entrySet()) {

            String baseType = entry.getKey(); // e.g. "N_2"
            List<Integer> atomIds = entry.getValue();

            // List of unique charges for this base type
            List<Double> uniqueCharges = new ArrayList<>();

            for (int atomId : atomIds) {
                double q = chargeMap.get(atomId);

                boolean found = false;
                for (double uq : uniqueCharges) {
                    if (Math.abs(q - uq) < Q_TOL) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    uniqueCharges.add(q);
                }
            }

            // Sort charges for deterministic q1, q2, ...
            uniqueCharges.sort(Double::compareTo);

            // Map charge → qIndex
            Map<Double, Integer> chargeIndexMap = new HashMap<>();
            for (int i = 0; i < uniqueCharges.size(); i++) {
                chargeIndexMap.put(uniqueCharges.get(i), i + 1);
            }

            for (int atomId : atomIds) {
                double q = chargeMap.get(atomId);

                int qIndex = -1;
                for (int i = 0; i < uniqueCharges.size(); i++) {
                    if (Math.abs(q - uniqueCharges.get(i)) < Q_TOL) {
                        qIndex = i + 1;
                        break;
                    }
                }

                if (qIndex < 0) {
                    throw new IllegalStateException(
                            "Charge " + q + " not matched for baseType " + baseType
                    );
                }

                String expandedType = baseType + "_q" + qIndex;
                expandedAtomTypeMap.put(atomId, expandedType);
            }

        }
    }

    public void advancedChargeMap(Map<Integer, String> expandedAtomTypeMap, Map<Integer, Double> chargeMap){
// Group atoms by base atom type
        Map<String, List<Integer>> atomsByBaseType = new HashMap<>();

        for (int atomId : atomTypeMap.keySet()) {
            String baseType = atomTypeMap.get(atomId);
            atomsByBaseType
                    .computeIfAbsent(baseType, k -> new ArrayList<>())
                    .add(atomId);
        }

// Process each base atom type
        for (Map.Entry<String, List<Integer>> entry : atomsByBaseType.entrySet()) {

            String baseType = entry.getKey();
            List<Integer> atomIds = entry.getValue();

            // Collect unique charges (with tolerance)
            List<Double> uniqueCharges = new ArrayList<>();

            for (int atomId : atomIds) {
                double q = chargeMap.get(atomId);

                boolean found = false;
                for (double uq : uniqueCharges) {
                    if (Math.abs(q - uq) < Q_TOL) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    uniqueCharges.add(q);
                }
            }

            // Sort for deterministic q1, q2, ...
            uniqueCharges.sort(Double::compareTo);

            // Assign q-indices and fill maps
            for (int atomId : atomIds) {

                double q = chargeMap.get(atomId);
                int qIndex = -1;
                double representativeQ = 0.0;

                for (int i = 0; i < uniqueCharges.size(); i++) {
                    double uq = uniqueCharges.get(i);
                    if (Math.abs(q - uq) < Q_TOL) {
                        qIndex = i + 1;
                        representativeQ = uq;
                        break;
                    }
                }

                String advancedType = baseType + "_q" + qIndex;
                expandedAtomTypeMap.put(atomId, advancedType);

                // Store charge once per advanced type
                advChargeMap.putIfAbsent(advancedType, representativeQ);
            }
        }
    }



    private void setUniqueElements(Set<String> uniqueElements){this.uniqueElements = uniqueElements;}

    public Map<Integer, String> getAtomMap(){return atomMap;}
    public Map<Integer, Vector> getPositionMap(){return positionMap;}
    public Set<String> getUniqueElements(){return uniqueElements;}
    /** Matches headers even if file has extra spaces like "Bond Coeffs " */
    private boolean equalsHeader(String line, String header) {
        return line.equalsIgnoreCase(header) || line.equalsIgnoreCase(header + ":");
    }

    /** If a line is mostly letters/spaces (not numeric data), treat as heading */
    private boolean looksLikeHeading(String line) {
        // headings like "Masses", "Pair Coeffs", "Bond Coeffs", etc.
        // numeric lines should start with digit or minus sign
        char c = line.charAt(0);
        if (Character.isDigit(c) || c == '-' || c == '+') return false;
        // If it contains at least one letter and not too many digits, call it a heading
        boolean hasLetter = false;
        for (int i = 0; i < line.length(); i++) {
            if (Character.isLetter(line.charAt(i))) { hasLetter = true; break; }
        }
        return hasLetter;
    }


    private enum Section {
        NONE, atoms,
        bond_Coeffs, angle_Coeffs, dihedral_Coeffs, torsion_Coeffs, improper_Coeffs,
        bonds, angles, torsions, dihedrals, impropers;
    }
    private void setConnectivity(ArrayList<ArrayList<Integer>>connectivityGrapehene){
        this.connectivityGrapehene = connectivityGrapehene;
    }
    private void setPotentials(Map<Integer, double[]> bondPotentials, Map<Integer, double[]> anglePotentials, Map<Integer, double[]> torsionPotentials,  Map<Integer, double[]> improperPotentials,  Map<Integer, double[]> dihedralPotentials, Map<Integer, Integer> elementNumMap, Map<Integer, Double> chargeMap, Map<Integer, String> atomMap, Map<Integer, Vector> positionMap){
        this.bondPotential = bondPotentials;
        this.anglePotential = anglePotentials;
        this.torsionPotential = torsionPotentials;
        this.atomMap = atomMap;
        this.positionMap = positionMap;
        this.chargeMap = chargeMap;
        this.elementNumMap = elementNumMap;
    }

    public void speciesGraphene (Map<Integer, Vector> listofPositions, Map<Integer, String> atomMap){
        ISpecies speciesGraphene ;
        for(int i = 0; i < atomMap.size(); i++){
            String atomName = atomMap.get(i);
        }
    }

    private String elementFromMass(double mass) {
        double tol = 1e-3;
        if (Math.abs(mass - 12.011) < tol) return "C";
        if (Math.abs(mass - 1.008)  < tol) return "H";
        if (Math.abs(mass - 15.999) < tol) return "O";
        if (Math.abs(mass - 14.007) < tol) return "N";
        if (Math.abs(mass - 40.078) < tol) return "Ca";
        if (Math.abs(mass - 30.974) < tol) return "P";
        throw new IllegalArgumentException("Unknown mass: " + mass);
    }
    private Map<Integer, String> makeUniqueAtomTypes(Map<Integer, Double> massByTypeId) {
        // -----------------------------------------
        // Step 1: find contiguous range starting at 1
        // -----------------------------------------
        List<Integer> typeIds = new ArrayList<>(massByTypeId.keySet());
        Collections.sort(typeIds);

        int maxContiguous = 0;
        for (int id : typeIds) {
            if (id == maxContiguous + 1) {
                maxContiguous = id;
            } else {
                break;
            }
        }
        // -----------------------------------------
        // Step 2: build unique atom types only for 1..maxContiguous
        // -----------------------------------------
        Map<String, Integer> elementCount = new HashMap<>();
        Map<Integer, String> out = new LinkedHashMap<>();

        for (int typeId = 1; typeId <= maxContiguous; typeId++) {
            Double mass = massByTypeId.get(typeId);
            if (mass == null) continue; // safety

            String element = elementFromMass(mass);

            int n = elementCount.getOrDefault(element, 0) + 1;
            elementCount.put(element, n);

            out.put(typeId, element + "_" + n);
        }
        return out;
    }
    private Map<Integer, String> makeAtomMap(Map<Integer, String > uniqueTypes, Map<Integer, Integer>elementNumMap){
        Map<Integer, String> atomMap = new HashMap<>();
        int num = 0;
        String atomName ="";
        for (int i = 1; i < elementNumMap.size()+1; i++){
            atomName = uniqueTypes.get(elementNumMap.get(i));
            atomMap.put(i, atomName);
        }
        return atomMap;
    }

    public ISpecies getSpecies(String confName, Vector vectorCOM, boolean isInfinite){
        AtomType typeNew;
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        GrapheneReaderXYZPDB grapheneReaderXYZPDB = new GrapheneReaderXYZPDB();
        AtomTypeFactory atomTypeFactory = new AtomTypeFactory();
        grapheneReaderXYZPDB.readDatafile(confName);
        positionsGraphene = grapheneReaderXYZPDB.getPositionMap();
        atomMapGraphene = grapheneReaderXYZPDB.getAtomMap();
        Space space = Space3D.getInstance();
        Vector center = space.makeVector();
        Vector dr = Vector.d(center.getD());

        for(int i = 1; i < atomMapGraphene.size()+1; i++) {
            String symbol = String.valueOf(atomMapGraphene.get(i));
            AtomType newName = atomTypeFactory.fromLabel(symbol);
            if (typeMapNew.containsKey(newName.getName())) {
                typeNew = typeMapNew.get(newName.getName());
            } else {
                typeNew = newName;
                typeMapNew.put(newName.getName(), typeNew);
            }

            Vector position = positionsGraphene.get(i);
            //System.out.println(position +  " " + typeNew + "  " + i);
            speciesBuilderNew.addAtom(typeNew, position,  "");
        }
        //System.out.println(typeMapNew + " typeMapNew");
        species= speciesBuilderNew.setDynamic(true).build();
        //   System.out.println(species.getMass() + " first");
        return species;
    }

    private void makeIntArrayList(GrapheneReaderXYZPDB grapheneReaderXYZPDB){
        IntArrayList[] hydroxylIntArr = new IntArrayList[50];
        IntArrayList[] carboxylIntArr = new IntArrayList[50];
        IntArrayList[] epoxyIntArr = new IntArrayList[50];
        Map<Integer, List<int[]>>mapTypeMap = grapheneReaderXYZPDB.getAngleTypeMap();
        for (int i = 0; i < mapTypeMap.size(); i ++){
            String combination = mapTypeMap.toString();
            String[] combinationSlices = combination.split("-");
            if (combinationSlices[0].equals("C") && combinationSlices[1].equals("O") && combinationSlices[2].equals("H")){
                List<int[]> angles = mapTypeMap.get(combination);
                hydroxylIntArr = grapheneReaderXYZPDB.makeIntArrFromList(angles);
            }else if(combinationSlices[0].equals("C") && combinationSlices[1].equals("C") && combinationSlices[2].equals("O")){
                List<int[]> angles = mapTypeMap.get(combination);
                carboxylIntArr = grapheneReaderXYZPDB.makeIntArrFromList(angles);
            }else if(combinationSlices[0].equals("C") && combinationSlices[1].equals("O") && combinationSlices[2].equals("C")){
                List<int[]> angles = mapTypeMap.get(combination);
                epoxyIntArr = grapheneReaderXYZPDB.makeIntArrFromList(angles);
            }
        }
        setIntArr(carboxylIntArr, hydroxylIntArr, epoxyIntArr);
    }

    private IntArrayList[] makeIntArrFromList(List<int[]> angles){
        IntArrayList[] result = new IntArrayList[angles.size()];
        int[] angle = new int[angles.get(0).length];
        for (int i = 0; i < angles.size(); i++) {
            angle = angles.get(i);
            IntArrayList ial = new IntArrayList(angle.length);
            for (int v : angle) {
                ial.add(v);
            }
            result[i] = ial;
        }
        return result;
    }

    public AtomType atomFactory(String atomType){
        AtomType atomType1 = new AtomType(Carbon.INSTANCE);
        return atomType1;
    }
    private void setBondTypeMap(Map<Integer, List<int[]>> bondTypeMap){this.bondTypeMap = bondTypeMap;}
    private void setAngleTypeMap(Map<Integer, List<int[]>> angleTypeMap){this.angleTypeMap = angleTypeMap;}
    private void setTorsionTypeMap(Map<Integer, List<int[]>> torsionTypeMap){this.torsionTypeMap = torsionTypeMap;}
    private void setImproperTypeMap(Map<Integer, List<int[]>> improperTypeMap){this.improperTypeMap = improperTypeMap;}
    private void setDihedralTypeMap(Map<Integer, List<int[]>> diheadralTypeMap){this.diheadralTypeMap = diheadralTypeMap;}
    public void setUniqueAtomTypeMap(Map<Integer, String> map) {this.uniqueAtomTypeMap = map;}
    private void setIntArr(IntArrayList[] carboxyIntArr, IntArrayList[] hydroxyIntArr, IntArrayList[] epoxyIntArr){
        this.epoxyIntArr = epoxyIntArr;
        this.hydroxyIntArr = hydroxyIntArr;
        this.carboxyIntArr = carboxyIntArr;
    }

    public Map<Integer, double[]> getBondPotential(){return bondPotential;}
    public Map<Integer, double[]> getAnglePotential(){return anglePotential;}
    public Map<Integer, double[]> getTorsionPotential(){return torsionPotential;}
    public Map<Integer, List<int[]>> getBondTypeMap(){return bondTypeMap;}
    public Map<Integer, List<int[]>> getAngleTypeMap(){return angleTypeMap;}
    public Map<Integer, List<int[]>> getTorsionTypeMap(){return torsionTypeMap;}
    public Map<Integer, Double> getChargeMap(){return chargeMap;}
    public Map<Integer, String> getUniqueAtomTypeMap() {return uniqueAtomTypeMap;}

    public static void main(String[] args) {
        GrapheneReaderXYZPDB grapheneReaderXYZPDB = new GrapheneReaderXYZPDB();
        grapheneReaderXYZPDB.readDatafile("D:\\Sem-X\\GO\\graphitis\\GO_sheet");
     //   grapheneReaderXYZPDB.readXYZfile("D:\\Sem-X\\GO\\graphitis\\GO_sheet");
        Map<Integer, String> atomMap = grapheneReaderXYZPDB.getAtomMap();
        Map<Integer, Vector> positions =grapheneReaderXYZPDB.getPositionMap();
        Map<Integer, String> uniqueAtomTypeMap = grapheneReaderXYZPDB.getUniqueAtomTypeMap();
        Map<Integer, Double> chargeMap = grapheneReaderXYZPDB.getChargeMap();
        System.out.println(atomMap);
        System.out.println(positions);
        System.out.println(uniqueAtomTypeMap);
        System.out.println(chargeMap);
       /* Map<Integer, double[]> bondPotentials = grapheneReaderXYZPDB.getBondPotential();
        Map<Integer, List<int[]>> bondTypeMap = grapheneReaderXYZPDB.getBondTypeMap();
        System.out.println(bondPotentials);
        System.out.println(bondTypeMap);*/
    }
}
