package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.potential.*;
import etomica.potential.OPLS_AA.BondTorsionOPLS;
import etomica.potential.UFF.P2HarmonicUFF;
import etomica.potential.UFF.P4BondInversionUFF;
import etomica.potential.UFF.P4BondTorsionUFF;
import etomica.potential.UFF.UFF;
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

public class GasOPLS {
    private ISpecies species;
    double massSum;
    public Map<Integer, String> atomMapGas = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectivityGrapehene = new ArrayList<>();
    public Map<Integer, String> modifiedAtomIdentifierMapGraphene = new HashMap<>();
    public Map<Integer, Vector> positionsGas = new HashMap<>();
    public Map<Integer,List<int[]>> bondTypeMap = new HashMap<>();
    public Map<Integer, double[]> pairCoeffsMap = new HashMap<>();
    public Map<Integer,List<int[]>> angleTypeMap = new HashMap<>();
    public Map<Integer,List<int[]>> torsionTypeMap = new HashMap<>();
    public Map<Integer,List<int[]>> diheadralTypeMap = new HashMap<>();
    public Map<Integer,List<int[]>> improperTypeMap = new HashMap<>();
    public Map<String, AtomType> typeMap = new HashMap<>();
    public List<int[]> dupletsSorted, tripletsSorted, quadrupletsSorted = new ArrayList<>();
    Map<String, AtomType> typeMapNew = new HashMap<>();
    public Map<Integer, String> atomMap = new HashMap<>();
    protected Set<String> uniqueElements = new HashSet<>();
    public Map<Integer, String> uniqueAtomTypeMap;
    public Map<Integer,double[]> bondPotential = new HashMap<>();
    public Map<Integer,double[]> anglePotential = new HashMap<>();
    public Map<Integer,double[]> torsionPotential = new HashMap<>();
    public Map<String, double[]> coeffPotential = new HashMap<>();
    public Map<Integer, Integer> elementNumMap = new HashMap<>();
    public Map<Integer, Double> chargeMap = new HashMap<>();
    public Map<Integer, Vector> positionMap = new HashMap<>();
    double Q_TOL = 0.0001;
    protected IntArrayList[] hydroxyIntArr, epoxyIntArr, carboxyIntArr;

    Map<Integer, String> expandedAtomTypeMap = new HashMap<>();
    Map<String, Double> chargeCoeff = new HashMap<>();
    public void readOPLSDataFile(String filePath){
        String fileName = filePath + ".lmp";   // <-- you said .data file
        Map<Integer, double[]> bondPotentials  = new HashMap<>();
        Map<Integer, double[]> anglePotentials = new HashMap<>();
        Map<Integer, double[]> torsionPotentials = new HashMap<>();
        Map<Integer, double[]> improperPotentials = new HashMap<>();
        Map<Integer, double[]> dihedralPotentials = new HashMap<>();
        Map<Integer,List<int[]>> improperTypeMap = new HashMap<>();
        Map<String, double[]> coeffPotential = new HashMap<>();
        Map<String, Double> chargePotential = new HashMap<>();
        boolean inMasses = false;
        Section section = Section.NONE;
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
                if (equalsHeader(line, "Pair Coeffs"))     { section = Section.pairCoeff;     continue; }
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
                    case pairCoeff: {
                        String[] t = line.split("\\s+");
                        double[] doubles = new double[2];
                        // atomID   moleculeID  atomTypeID  partialCharge   X   Y   Z
                        if (t.length >= 2){
                            int numZero = Integer.parseInt(t[0]);
                            double numTwo = Double.parseDouble(t[1]);
                            double numThree = Double.parseDouble(t[2]);
                            doubles = new double[]{numTwo, numThree};
                            pairCoeffsMap.put(numZero, doubles);
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

            GrapheneReaderXYZPDB grapheneReaderXYZPDB = new GrapheneReaderXYZPDB();
            setPositions(positionMap);
            setBondTypeMap(bondTypeMap);
            setAngleTypeMap(angleTypeMap);
            setTorsionTypeMap(torsionTypeMap);
            setDihedralTypeMap(diheadralTypeMap);
            setImproperTypeMap(improperTypeMap);
            setConnectivity(connectivityGrapehene);
            Map<Integer, String> uniqueTypes = grapheneReaderXYZPDB.makeUniqueAtomTypes(massByTypeId);
            setUniqueAtomTypeMap(uniqueTypes);
            coeffPotential = makeCoeffPotential(uniqueTypes, pairCoeffsMap, chargeMap);
            atomMap = grapheneReaderXYZPDB.makeAtomMap(uniqueTypes, elementNumMap);
            setPotentials(bondPotentials, anglePotentials, torsionPotentials, improperPotentials, dihedralPotentials, coeffPotential, elementNumMap, chargeMap, atomMap, positionMap);
            System.out.println("Done Reading: " + fileName);

        } catch (IOException e) {
            throw new RuntimeException("Problem reading from " + fileName + ", caught IOException: " + e.getMessage(), e);
        }
    }
    private enum Section {
        NONE, atoms, pairCoeff,
        bond_Coeffs, angle_Coeffs, dihedral_Coeffs, torsion_Coeffs, improper_Coeffs,
        bonds, angles, torsions, dihedrals, impropers;
    }

    public Map<String, double[]> makeCoeffPotential(Map<Integer, String> uniqueAtomTypeMap, Map<Integer, double[]> coeffPairs, Map<Integer, Double> chargeMap){
        Map<String, double[]> coeffPotentialMap = new HashMap<>();
        Map<String, Double> chargeCoeff = new HashMap<>();
        double[] doubles = new double[2];
        double chargeVal;
        String atomType = "";
        for (int i=1; i < uniqueAtomTypeMap.size()+1; i++ ){
            doubles = coeffPairs.get(i);
            atomType = uniqueAtomTypeMap.get(i);
            coeffPotentialMap.put(atomType, doubles);

            chargeVal = chargeMap.get(i);
            chargeCoeff.put(uniqueAtomTypeMap.get(i), chargeVal);
        }
        setChargeCoeff(chargeCoeff);
        return coeffPotentialMap;
    }
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

    public ISpecies getSpecies(String confName, Vector vectorCOM, boolean isInfinite){
        AtomType typeNew;
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        GasOPLS gasOPLS = new GasOPLS();
        AtomTypeFactory atomTypeFactory = new AtomTypeFactory();
        gasOPLS.readOPLSDataFile(confName);
        positionMap = gasOPLS.getPosition();
        atomMapGas = gasOPLS.getAtomMap();
        coeffPotential = gasOPLS.getCoeffPotential();
        chargeCoeff = gasOPLS.getChargeCoeff();
        Space space = Space3D.getInstance();
        Vector center = space.makeVector();
        Vector dr = Vector.d(center.getD());

        for(int i = 1; i < atomMapGas.size()+1; i++) {
            String symbol = String.valueOf(atomMapGas.get(i));
            AtomType newName = atomTypeFactory.fromLabel(symbol);
            if (typeMapNew.containsKey(newName.getName())) {
                typeNew = typeMapNew.get(newName.getName());
            } else {
                typeNew = newName;
                typeMapNew.put(newName.getName(), typeNew);
            }

            Vector position = positionMap.get(i);
            //System.out.println(position +  " " + typeNew + "  " + i);
            speciesBuilderNew.addAtom(typeNew, position,  "");
        }
        //System.out.println(typeMapNew + " typeMapNew");
        species= speciesBuilderNew.setDynamic(true).build();
        //   System.out.println(species.getMass() + " first");
        return species;
    }

    public void makePotential(String gasName, PotentialMasterBonding pmBonding){
        GasOPLS gasOPLS = new GasOPLS();

        Map<Integer, List<int[]>> bondTypesMap2 = new HashMap<>();
        Map<Integer, List<int[]>> angleTypesMap2 = new HashMap<>();
        Map<Integer, List<int[]>> torsionTypesMap2 = new HashMap<>();

        ISpecies species1 = gasOPLS.getSpecies(gasName, new Vector3D(0,0,0),false);
        bondTypesMap2 = gasOPLS.getBondTypeMap();
        angleTypesMap2 = gasOPLS.getAngleTypeMap();
        torsionTypesMap2 = gasOPLS.getTorsionTypeMap();
        Map<Integer, double[]> bondPotential = gasOPLS.getBondPotential();
        Map<Integer, double[]> anglePotential = gasOPLS.getBondPotential();
        Map<Integer, double[]> torsionPotential = gasOPLS.getBondPotential();

        for (int i = 1; i < bondTypesMap2.size(); i++){
            double[] potential = bondPotential.get(i);
            P2HarmonicUFF p2HarmonicUFF = new P2HarmonicUFF(potential[0], potential[1]);
            pmBonding.setBondingPotentialPair(species1, p2HarmonicUFF, bondTypesMap2.get(i));
        }

        for (int i = 1; i < angleTypesMap2.size(); i++){
            double[] potential = anglePotential.get(i);
            P3BondAngle p3BondAngle = new P3BondAngle(potential[0], potential[1] );
            pmBonding.setBondingPotentialTriplet(species1, p3BondAngle, angleTypesMap2.get(i));
        }

        for (int i =1; i < torsionTypesMap2.size(); i++){
            double[] potential = torsionPotential.get(i);
            BondTorsionOPLS bondTorsionOPLS = new BondTorsionOPLS(potential[0], potential[1], potential[2], potential[3]);
            pmBonding.setBondingPotentialQuad(species1, bondTorsionOPLS, angleTypesMap2.get(i));
        }
    }

    private void setPositions(Map<Integer, Vector> positions){this.positionMap = positions;}

    private boolean equalsHeader(String line, String header) {
        return line.equalsIgnoreCase(header) || line.equalsIgnoreCase(header + ":");
    }
    private void setBondTypeMap(Map<Integer, List<int[]>> bondTypeMap){this.bondTypeMap = bondTypeMap;}
    private void setAngleTypeMap(Map<Integer, List<int[]>> angleTypeMap){this.angleTypeMap = angleTypeMap;}
    private void setTorsionTypeMap(Map<Integer, List<int[]>> torsionTypeMap){this.torsionTypeMap = torsionTypeMap;}
    private void setImproperTypeMap(Map<Integer, List<int[]>> improperTypeMap){this.improperTypeMap = improperTypeMap;}
    private void setDihedralTypeMap(Map<Integer, List<int[]>> diheadralTypeMap){this.diheadralTypeMap = diheadralTypeMap;}
    public void setUniqueAtomTypeMap(Map<Integer, String> map) {this.uniqueAtomTypeMap = map;}
    public void setChargeCoeff(Map<String, Double> chargeCoeff){this.chargeCoeff = chargeCoeff;}

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
    public Map<String, Double> getChargeCoeff(){return chargeCoeff;}
    public Map<Integer, String> getUniqueAtomTypeMap() {return uniqueAtomTypeMap;}
    public Map<String, double[]> getCoeffPotential(){return coeffPotential;}
    public Map<Integer, Vector> getPosition(){return positionMap;}
    public Map<Integer, String> getAtomMap(){return atomMap;}
    private void setConnectivity(ArrayList<ArrayList<Integer>>connectivityGrapehene){
        this.connectivityGrapehene = connectivityGrapehene;
    }
    private void setPotentials(Map<Integer, double[]> bondPotentials, Map<Integer, double[]> anglePotentials, Map<Integer, double[]> torsionPotentials,  Map<Integer, double[]> improperPotentials,  Map<Integer, double[]> dihedralPotentials,Map<String, double[]> coeffPotential, Map<Integer, Integer> elementNumMap, Map<Integer, Double> chargeMap, Map<Integer, String> atomMap, Map<Integer, Vector> positionMap){
        this.bondPotential = bondPotentials;
        this.anglePotential = anglePotentials;
        this.torsionPotential = torsionPotentials;
        this.coeffPotential = coeffPotential;
        this.atomMap = atomMap;
        this.positionMap = positionMap;
        this.chargeMap = chargeMap;
        this.elementNumMap = elementNumMap;
    }

    public static void main(String[] args) {
        GasOPLS gasOPLS = new GasOPLS();
        gasOPLS.getSpecies("D:\\Sem-X\\GO\\gasOPLS\\1butene", new Vector3D(0,0,0),false);
    }
}
