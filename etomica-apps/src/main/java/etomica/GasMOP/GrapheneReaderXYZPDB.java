package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.potential.GAFF.P4BondTorsion1;
import etomica.potential.GAFF.P4BondTorsion2;
import etomica.potential.GAFF.P4BondTorsion3;
import etomica.potential.OPLS_AA.BondTorsionOPLS;
import etomica.potential.UFF.*;
import etomica.potential.compute.PotentialComputeEwaldFourier;
import etomica.potential.ewald.P2Ewald1Real;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.units.*;
import etomica.util.collections.IntArrayList;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class GrapheneReaderXYZPDB {
    private ISpecies species;
    double massSum;
    public Map<Integer, String> atomMapGraphene = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectivityGrapehene = new ArrayList<>();
    public Map<Integer, String> modifiedAtomIdentifierMapGraphene = new HashMap<>();
    public Map<Integer, Vector> positionsGraphene = new HashMap<>();
    public Map<Integer,List<int[]>> bondTypeMap = new HashMap<>();
    public Map<Integer,List<int[]>> angleTypeMap = new HashMap<>();
    public Map<Integer, double[]> pairCoeffsMap = new HashMap<>();
    public Map<Integer,List<int[]>> torsionTypeMap = new HashMap<>();
    public Map<Integer,List<int[]>> improperTypeMap = new HashMap<>();
    public Map<String, AtomType> typeMap = new HashMap<>();
    Map<String, AtomType> typeMapNew = new HashMap<>();
    public Map<Integer, String> atomMap = new HashMap<>();
    public Set<String> uniqueElements = new HashSet<>();
    public Map<Integer, String> uniqueAtomTypeMap;
    public Map<String, Double> chargeCoeff = new HashMap<>();
    public Map<String, double[]> coeffPotential = new HashMap<>();
    public Map<Integer,double[]> bondPotential = new HashMap<>();
    public Map<Integer,double[]> anglePotential = new HashMap<>();
    public Map<Integer,double[]> torsionPotential = new HashMap<>();
    public Map<Integer,double[]> improperPotential = new HashMap<>();
    public Map<Integer, Integer> elementNumMap = new HashMap<>();
    public  Map<String, double[]> coeffPotentialMapNew = new HashMap<>();
    public Map<String, List<Integer>> newAtomTypeMap = new HashMap<>();
    public Map<String, Double> newTypeChargeMap = new HashMap<>();
    public Map<String, Map<Double, String>> typeChargeMap = new HashMap<>();
    public Map<Integer, String> newAtomMap = new HashMap<>();
    public Map<Integer, Double> chargeMap = new HashMap<>();
    public Map<Integer, Vector> positionMap = new HashMap<>();
    double Q_TOL = 0.0001;
    protected IntArrayList[] hydroxyIntArr, epoxyIntArr, carboxyIntArr;

    public Map<Integer, String> atomTypeMap;
    public Map<Integer, String> expandedAtomTypeMap = new HashMap<>();
    public Map<String, Double> advChargeMap = new HashMap<>();

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
        GrapheneReaderXYZPDB grapheneReaderXYZPDB = new GrapheneReaderXYZPDB();
        // ---- choose correct extension ----
        String fileName = confName + ".data";   // <-- you said .data file (not .xyz)
        // If you really want: confName+".xyz" then change back.

        // ---- your maps (use your existing ones if already fields) ----
        // If these are fields already, DELETE these local declarations.
        Map<Integer, double[]> bondPotentials  = new HashMap<>();
        Map<Integer, double[]> anglePotentials = new HashMap<>();
        Map<Integer, double[]> torsionPotentials = new HashMap<>();
        Map<Integer, double[]> improperPotentials = new HashMap<>();
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
                if (equalsHeader(line, "Pair Coeffs"))     { section = Section.pairCoeff;    continue; }
                if (equalsHeader(line, "Bond Coeffs"))     { section = Section.bond_Coeffs;    continue; }
                if (equalsHeader(line, "Angle Coeffs"))    { section = Section.angle_Coeffs;   continue; }
                if (equalsHeader(line, "Dihedral Coeffs")) { section = Section.torsion_Coeffs;continue; }
                if (equalsHeader(line, "Improper Coeffs")) { section = Section.improper_Coeffs;continue; }

                // Topology sections (optional if your .data has them)
                if (equalsHeader(line, "Atoms"))     { section = Section.atoms;     continue; }
                if (equalsHeader(line, "Bonds"))     { section = Section.bonds;     continue; }
                if (equalsHeader(line, "Angles"))    { section = Section.angles;    continue; }
                if (equalsHeader(line, "Torsions"))  { section = Section.torsions;  continue; }
                if (equalsHeader(line, "Dihedrals")) { section = Section.torsions; continue; }
                if (equalsHeader(line, "Impropers")) { section = Section.impropers; continue; }

                // If line is another random heading (letters but not numeric), exit section safely
                if (looksLikeHeading(line)) {
                    section = Section.NONE;
                    continue;
                }

                // -------------------------
                // parse based on section
                // -------------------------
           //     System.out.println(line);
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

                    case improper_Coeffs: {
                        String[] t = line.split("\\s+");
                        if (t.length >= 1) {
                            int id = Integer.parseInt(t[0]);
                            double[] params = new double[t.length - 1];
                            for (int i = 1; i < t.length; i++) {
                                params[i - 1] = Double.parseDouble(t[i]);
                            }
                            improperPotentials.put(id, params);
                        }
                        break;
                    }

                    case torsion_Coeffs: {
                        String[] t = line.split("\\s+");
                        if (t.length >= 3) {
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
                            int atom1    = Integer.parseInt(t[2])-1;
                            int atom2    = Integer.parseInt(t[3])-1;

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
                            int atom1    = Integer.parseInt(t[2])-1;
                            int atom2    = Integer.parseInt(t[3])-1;
                            int atom3    = Integer.parseInt(t[4])-1;

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
                            int atom1    = Integer.parseInt(t[2])-1;
                            int atom2    = Integer.parseInt(t[3])-1;
                            int atom3    = Integer.parseInt(t[4])-1;
                            int atom4    = Integer.parseInt(t[5])-1;
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
            GasOPLS gasOPLS = new GasOPLS();
            setBondTypeMap(bondTypeMap);
            setAngleTypeMap(angleTypeMap);
            setTorsionTypeMap(torsionTypeMap);
            setImproperTypeMap(improperTypeMap);
            setConnectivity(connectivityGrapehene);
            Map<Integer, String> uniqueTypes = makeUniqueAtomTypes(massByTypeId);
            Map<Integer, String> usedUniqueTypes = new LinkedHashMap<>();
            for (Integer typeId : new TreeSet<>(elementNumMap.values())) {   // sorted unique ids
                String name = uniqueTypes.get(typeId);
                if (name != null) usedUniqueTypes.put(typeId, name);
            }
            System.out.println("usedUniqueTypes " + usedUniqueTypes);
           // Map<String, double[]> coeffPotential1 = gasOPLS.makeCoeffPotential(uniqueTypes, pairCoeffsMap, chargeMap);
         //   System.out.println("pot1 "+coeffPotential1);

            Map<Integer, String> newAtomMap = new HashMap<>();
            Map<String, double[]> coeffPotentialMapNew = new HashMap<>();
            Map<String, List<Integer>> newAtomTypeMap = new HashMap<>();
            Map<String, Double> newTypeChargeMap = new HashMap<>();
            Map<String, Map<Double, String>> typeChargeMap = new HashMap<>();
            setUniqueAtomTypeMap(usedUniqueTypes);
            atomMap = makeAtomMap(uniqueTypes, elementNumMap);
            grapheneReaderXYZPDB.makeCoeffPotential(usedUniqueTypes, pairCoeffsMap, chargeMap, atomMap, newAtomMap, coeffPotentialMapNew, newAtomTypeMap, newTypeChargeMap, typeChargeMap);
            setNewMaps(newAtomMap, coeffPotentialMapNew, newAtomTypeMap, newTypeChargeMap, typeChargeMap);
            setPotentials(bondPotentials, anglePotentials, torsionPotentials, improperPotentials, elementNumMap, positionMap);
            System.out.println("Done Reading: " + fileName);

        } catch (IOException e) {
            throw new RuntimeException("Problem reading from " + fileName + ", caught IOException: " + e.getMessage(), e);
        }
    }


    public void readDatafileData(String confName) {
        GrapheneReaderXYZPDB grapheneReaderXYZPDB = new GrapheneReaderXYZPDB();
        // ---- choose correct extension ----
        String fileName = confName + ".data";   // <-- you said .data file (not .xyz)
        // If you really want: confName+".xyz" then change back.

        // ---- your maps (use your existing ones if already fields) ----
        // If these are fields already, DELETE these local declarations.
        Map<Integer, double[]> bondPotentials  = new HashMap<>();
        Map<Integer, double[]> anglePotentials = new HashMap<>();
        Map<Integer, double[]> torsionPotentials = new HashMap<>();
        Map<Integer, double[]> improperPotentials = new HashMap<>();
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
                if (equalsHeader(line, "Pair Coeffs"))     { section = Section.pairCoeff;    continue; }
                if (equalsHeader(line, "Bond Coeffs"))     { section = Section.bond_Coeffs;    continue; }
                if (equalsHeader(line, "Angle Coeffs"))    { section = Section.angle_Coeffs;   continue; }
                if (equalsHeader(line, "Dihedral Coeffs")) { section = Section.torsion_Coeffs;continue; }
                if (equalsHeader(line, "Improper Coeffs")) { section = Section.improper_Coeffs;continue; }

                // Topology sections (optional if your .data has them)
                if (equalsHeader(line, "Atoms"))     { section = Section.atoms;     continue; }
                if (equalsHeader(line, "Bonds"))     { section = Section.bonds;     continue; }
                if (equalsHeader(line, "Angles"))    { section = Section.angles;    continue; }
                if (equalsHeader(line, "Torsions"))  { section = Section.torsions;  continue; }
                if (equalsHeader(line, "Dihedrals")) { section = Section.torsions; continue; }
                if (equalsHeader(line, "Impropers")) { section = Section.impropers; continue; }

                // If line is another random heading (letters but not numeric), exit section safely
                if (looksLikeHeadingAmber(line)) {
                    section = Section.NONE;
                    continue;
                }

                // -------------------------
                // parse based on section
                // -------------------------

                switch (section) {

                    case bond_Coeffs: {
                        String clean = line.split("#", 2)[0].trim();
                        if (clean.isEmpty()) break;

                        String[] t = clean.split("\\s+");

                        if (t.length >= 4 && t[0].equalsIgnoreCase("bond_coeff")) {
                            int id = Integer.parseInt(t[1]);
                            double k  = Double.parseDouble(t[2]);
                            double r0 = Double.parseDouble(t[3]);
                            bondPotentials.put(id, new double[]{k, r0});
                        } else if (t.length >= 3 && Character.isDigit(t[0].charAt(0))) {
                            // fallback for plain numeric style:  id k r0
                            int id = Integer.parseInt(t[0]);
                            double k  = Double.parseDouble(t[1]);
                            double r0 = Double.parseDouble(t[2]);
                            bondPotentials.put(id, new double[]{k, r0});
                        } else {
                            throw new RuntimeException("Bad bond coeff line: " + Arrays.toString(t));
                        }
                        break;
                    }

                    case angle_Coeffs: {
                        String clean = line.split("#", 2)[0].trim();
                        if (clean.isEmpty()) break;

                        String[] t = clean.split("\\s+");

                        if (t.length >= 4 && t[0].equalsIgnoreCase("angle_coeff")) {
                            int id = Integer.parseInt(t[1]);
                            double k = Double.parseDouble(t[2]);
                            double theta0 = Double.parseDouble(t[3]);
                            anglePotentials.put(id, new double[]{k, theta0});
                        } else if (t.length >= 3 && Character.isDigit(t[0].charAt(0))) {
                            // fallback for plain numeric style: id k theta0
                            int id = Integer.parseInt(t[0]);
                            double k = Double.parseDouble(t[1]);
                            double theta0 = Double.parseDouble(t[2]);
                            anglePotentials.put(id, new double[]{k, theta0});
                        } else {
                            throw new RuntimeException("Bad angle coeff line: " + Arrays.toString(t));
                        }
                        break;
                    }

                    case pairCoeff: {
                        String clean = line.split("#", 2)[0].trim();
                        if (clean.isEmpty()) break;

                        String[] t = clean.split("\\s+");

                        if (t.length >= 5 && t[0].equalsIgnoreCase("pair_coeff")) {
                            int type1 = Integer.parseInt(t[1]);
                            int type2 = Integer.parseInt(t[2]);
                            double epsilon = Double.parseDouble(t[3]);
                            double sigma = Double.parseDouble(t[4]);

                            // store only self coefficients like pair_coeff 8 8 ...
                            if (type1 == type2) {
                                pairCoeffsMap.put(type1, new double[]{epsilon, sigma});
                            }
                        } else {
                            throw new RuntimeException("Bad pair coeff line: " + Arrays.toString(t));
                        }
                        break;
                    }

                    case improper_Coeffs: {
                        String clean = line.split("#", 2)[0].trim();
                        if (clean.isEmpty()) break;

                        String[] t = clean.split("\\s+");

                        if (t.length >= 2 && t[0].equalsIgnoreCase("improper_coeff")) {
                            int id = Integer.parseInt(t[1]);
                            double[] params = new double[t.length - 2];
                            for (int i = 2; i < t.length; i++) {
                                params[i - 2] = Double.parseDouble(t[i]);
                            }
                            improperPotentials.put(id, params);
                        } else if (t.length >= 2 && Character.isDigit(t[0].charAt(0))) {
                            // fallback for plain numeric style: id p1 p2 ...
                            int id = Integer.parseInt(t[0]);
                            double[] params = new double[t.length - 1];
                            for (int i = 1; i < t.length; i++) {
                                params[i - 1] = Double.parseDouble(t[i]);
                            }
                            improperPotentials.put(id, params);
                        } else {
                            throw new RuntimeException("Bad improper coeff line: " + Arrays.toString(t));
                        }
                        break;
                    }

                    case torsion_Coeffs: {
                        String clean = line.split("#", 2)[0].trim();
                        if (clean.isEmpty()) break;

                        String[] t = clean.split("\\s+");

                        if (t.length >= 2 && t[0].equalsIgnoreCase("dihedral_coeff")) {
                            int id = Integer.parseInt(t[1]);
                            double[] params = new double[t.length - 2];
                            for (int i = 2; i < t.length; i++) {
                                params[i - 2] = Double.parseDouble(t[i]);
                            }
                            torsionPotentials.put(id, params);
                        } else if (t.length >= 2 && Character.isDigit(t[0].charAt(0))) {
                            // fallback for plain numeric style: id p1 p2 ...
                            int id = Integer.parseInt(t[0]);
                            double[] params = new double[t.length - 1];
                            for (int i = 1; i < t.length; i++) {
                                params[i - 1] = Double.parseDouble(t[i]);
                            }
                            torsionPotentials.put(id, params);
                        } else {
                            throw new RuntimeException("Bad dihedral coeff line: " + Arrays.toString(t));
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
                          //  chargeMap.put(numZero, numThree);
                            chargeMap.put(numZero, 0.0);
                        }else {
                            throw new RuntimeException("t "+ Arrays.toString(t));
                        }
                        break;
                    }

                    case bonds: {
                        String[] t = line.split("\\s+");
                        if (t.length >= 3) {
                            int bondType = Integer.parseInt(t[1]);
                            int atom1    = Integer.parseInt(t[2])-1;
                            int atom2    = Integer.parseInt(t[3])-1;

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
                            int atom1    = Integer.parseInt(t[2])-1;
                            int atom2    = Integer.parseInt(t[3])-1;
                            int atom3    = Integer.parseInt(t[4])-1;

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
                            int atom1    = Integer.parseInt(t[2])-1;
                            int atom2    = Integer.parseInt(t[3])-1;
                            int atom3    = Integer.parseInt(t[4])-1;
                            int atom4    = Integer.parseInt(t[5])-1;
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
            GasOPLS gasOPLS = new GasOPLS();
            setBondTypeMap(bondTypeMap);
            setAngleTypeMap(angleTypeMap);
            setTorsionTypeMap(torsionTypeMap);
            setImproperTypeMap(improperTypeMap);
            setConnectivity(connectivityGrapehene);
            Map<Integer, String> uniqueTypes = makeUniqueAtomTypes(massByTypeId);
            Map<Integer, String> usedUniqueTypes = new LinkedHashMap<>();
            for (Integer typeId : new TreeSet<>(elementNumMap.values())) {   // sorted unique ids
                String name = uniqueTypes.get(typeId);
                if (name != null) usedUniqueTypes.put(typeId, name);
            }
            System.out.println("usedUniqueTypes " + usedUniqueTypes);
            // Map<String, double[]> coeffPotential1 = gasOPLS.makeCoeffPotential(uniqueTypes, pairCoeffsMap, chargeMap);
            //   System.out.println("pot1 "+coeffPotential1);
           // setUniqueAtomTypeMap(usedUniqueTypes);

            Map<Integer, String> newAtomMap = new HashMap<>();
            Map<String, double[]> coeffPotentialMapNew = new HashMap<>();
            Map<String, List<Integer>> newAtomTypeMap = new HashMap<>();
            Map<String, Double> newTypeChargeMap = new HashMap<>();
            Map<String, Map<Double, String>> typeChargeMap = new HashMap<>();
            setUniqueAtomTypeMap(usedUniqueTypes);
            atomMap = makeAtomMap(uniqueTypes, elementNumMap);
            grapheneReaderXYZPDB.makeCoeffPotential(usedUniqueTypes, pairCoeffsMap, chargeMap, atomMap, newAtomMap, coeffPotentialMapNew, newAtomTypeMap, newTypeChargeMap, typeChargeMap);
            setNewMaps(newAtomMap, coeffPotentialMapNew, newAtomTypeMap, newTypeChargeMap, typeChargeMap);
            setPotentials(bondPotentials, anglePotentials, torsionPotentials, improperPotentials,  elementNumMap, positionMap);
            System.out.println("Done Reading: " + fileName);

        } catch (IOException e) {
            throw new RuntimeException("Problem reading from " + fileName + ", caught IOException: " + e.getMessage(), e);
        }
    }
    public void makeCoeffPotential(Map<Integer, String> uniqueAtomTypeMap, Map<Integer, double[]> coeffPairs, Map<Integer, Double> chargeMap, Map<Integer, String> atomMap, Map<Integer, String> newAtomMap,    Map<String, double[]> coeffPotentialMapNew,  Map<String, List<Integer>> newAtomTypeMap,  Map<String, Double> newTypeChargeMap, Map<String, Map<Double, String>> typeChargeMap ){
        Map<String, double[]> coeffPotentialMap = new HashMap<>();
        Map<String, Double> chargeCoeff = new HashMap<>();
        double[] doubles = new double[2];
        double chargeVal;
        String atomType = "";

        Map<String, Set<Double>> temp = new HashMap<>();

        for (Integer atomId : atomMap.keySet()) {

            String atomType2 = atomMap.get(atomId);
            Double charge = chargeMap.get(atomId);

            if (charge == null) continue;

            temp.computeIfAbsent(atomType2, k -> new HashSet<>()).add(charge);
        }

        for (Map.Entry<String, Set<Double>> entry : temp.entrySet()) {
            String atomType2 = entry.getKey();
            Set<Double> charges = entry.getValue();

            Map<Double, String> chargeIndex = new HashMap<>();

            int i = 0;
            for (Double charge : charges) {
                chargeIndex.put(charge, atomType2 + "_" + i);
                i++;
            }

            typeChargeMap.put(atomType2, chargeIndex);
        }

       // System.out.println(typeChargeMap);

        for (Map.Entry<String, Map<Double, String>> entry : typeChargeMap.entrySet()) {

            Map<Double, String> chargeToType = entry.getValue();

            for (Map.Entry<Double, String> e : chargeToType.entrySet()) {

                Double charge = e.getKey();
                String newType = e.getValue();

                newTypeChargeMap.put(newType, charge);
            }
        }

        for (Integer atomId : atomMap.keySet()) {

            atomType = atomMap.get(atomId);
            Double charge = chargeMap.get(atomId);

            Map<Double, String> chargeToType = typeChargeMap.get(atomType);
            if (chargeToType == null) continue;

            String newType = chargeToType.get(charge);
            if (newType == null) continue;

            newAtomTypeMap
                    .computeIfAbsent(newType, k -> new ArrayList<>())
                    .add(atomId);
        }

        for (Map.Entry<String, List<Integer>> entry : newAtomTypeMap.entrySet()) {

            String newType = entry.getKey();
            List<Integer> atoms = entry.getValue();

            for (Integer atomId : atoms) {
                newAtomMap.put(atomId, newType);
            }
        }

        for (Integer typeId : uniqueAtomTypeMap.keySet()) {
          //  System.out.println(typeId);
            doubles = coeffPairs.get(typeId);
            atomType = uniqueAtomTypeMap.get(typeId);
            coeffPotentialMap.put(atomType, doubles);
        }

        for (String newType : newAtomTypeMap.keySet()) {

            // remove the last "_index"
            int pos = newType.lastIndexOf('_');
            String baseType = newType.substring(0, pos);

            double[] params = coeffPotentialMap.get(baseType);

            if (params == null) {
                throw new RuntimeException("No pair coeffs for " + baseType);
            }

            coeffPotentialMapNew.put(newType, params);
        }
        //System.out.println(coeffPotentialMap);


    /*    for (String newType : newAtomTypeMap.keySet()) {

            // extract base atom type
            int pos = newType.lastIndexOf('_');
            String baseType = newType.substring(0, pos);

            // get LJ params of base type
            double[] params = coeffPotentialMap.get(baseType);

            System.out.println("------------------------------------------------");
            System.out.println("New atom type : " + newType);
            System.out.println("Base atom type: " + baseType);

            if (params == null) {
                System.out.println("ERROR: No coefficients found for " + baseType);
                continue;
            }

            System.out.println("Base coeffs   : " + Arrays.toString(params));

            // assign to new map
            coeffPotentialMapNew.put(newType, params.clone());

            System.out.println("Assigned coeff: " + Arrays.toString(coeffPotentialMapNew.get(newType)));
        }*/
       // setChargeCoeff(chargeCoeff);

        //return coeffPotentialMap;
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

    public Map<Integer, String> getAtomMap(){return newAtomMap;}
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

    private boolean looksLikeHeadingAmber(String line) {
        line = line.trim();
        if (line.isEmpty()) return false;

        // comment-only lines are not headings
        if (line.startsWith("#")) return false;

        // coefficient command lines inside sections are data, not headings
        String lower = line.toLowerCase();
        if (lower.startsWith("pair_coeff")) return false;
        if (lower.startsWith("bond_coeff")) return false;
        if (lower.startsWith("angle_coeff")) return false;
        if (lower.startsWith("dihedral_coeff")) return false;
        if (lower.startsWith("improper_coeff")) return false;

        // ordinary numeric data lines are not headings
        char c = line.charAt(0);
        if (Character.isDigit(c) || c == '-' || c == '+') return false;

        // recognized section headers
        if (equalsHeader(line, "Masses")) return true;
        if (equalsHeader(line, "Pair Coeffs")) return true;
        if (equalsHeader(line, "Bond Coeffs")) return true;
        if (equalsHeader(line, "Angle Coeffs")) return true;
        if (equalsHeader(line, "Dihedral Coeffs")) return true;
        if (equalsHeader(line, "Improper Coeffs")) return true;
        if (equalsHeader(line, "Atoms")) return true;
        if (equalsHeader(line, "Bonds")) return true;
        if (equalsHeader(line, "Angles")) return true;
        if (equalsHeader(line, "Dihedrals")) return true;
        if (equalsHeader(line, "Impropers")) return true;
        if (equalsHeader(line, "Velocities")) return true;

        return false;
    }


    private enum Section {
        NONE, atoms,pairCoeff,
        bond_Coeffs, angle_Coeffs,torsion_Coeffs, improper_Coeffs,
        bonds, angles, torsions, impropers;
    }
    private void setConnectivity(ArrayList<ArrayList<Integer>>connectivityGrapehene){
        this.connectivityGrapehene = connectivityGrapehene;
    }

    public void setNewMaps(Map<Integer, String> newAtomMap, Map<String, double[]> coeffPotentialMapNew, Map<String, List<Integer>> newAtomTypeMap, Map<String, Double> newTypeChargeMap, Map<String, Map<Double, String>> typeChargeMap){
       this.newAtomMap = newAtomMap;
       this.coeffPotentialMapNew = coeffPotentialMapNew;
       this.newAtomTypeMap = newAtomTypeMap;
       this.newTypeChargeMap = newTypeChargeMap;
       this.typeChargeMap = typeChargeMap;
    }
    public void setPotentials(Map<Integer, double[]> bondPotentials, Map<Integer, double[]> anglePotentials, Map<Integer, double[]> torsionPotentials,  Map<Integer, double[]> improperPotentials,  Map<Integer, Integer> elementNumMap, Map<Integer, Vector> positionMap){
        this.bondPotential = bondPotentials;
        this.anglePotential = anglePotentials;
        this.torsionPotential = torsionPotentials;
        this.improperPotential = improperPotentials;
        this.positionMap = positionMap;
        this.elementNumMap = elementNumMap;
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
    public Map<Integer, String> makeUniqueAtomTypes(Map<Integer, Double> massByTypeId) {
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
    public Map<Integer, String> makeAtomMap(Map<Integer, String> uniqueTypes,
                                            Map<Integer, Integer> elementNumMap) {
        Map<Integer, String> atomMap = new HashMap<>();

        for (Map.Entry<Integer, Integer> e : elementNumMap.entrySet()) {
            int atomId = e.getKey();      // the atom index
            int typeId = e.getValue();    // the typeId (1..13)

            String atomName = uniqueTypes.get(typeId);
            if (atomName == null) {
                throw new IllegalArgumentException("Missing typeId " + typeId + " for atomId " + atomId);
            }
            atomMap.put(atomId, atomName);
        }
        return atomMap;
    }

    public ISpecies getSpecies(String confName, Vector vectorCOM, boolean isInfinite, boolean isAmber){
        AtomType typeNew;
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        GrapheneReaderXYZPDB grapheneReaderXYZPDB = new GrapheneReaderXYZPDB();
        AtomTypeFactory atomTypeFactory = new AtomTypeFactory();
        if (!isAmber){
            grapheneReaderXYZPDB.readDatafile(confName);
        }else {
            grapheneReaderXYZPDB.readDatafileData(confName);
        }

        positionsGraphene = grapheneReaderXYZPDB.getPositionMap();
        atomMapGraphene = grapheneReaderXYZPDB.getAtomMap();
        newTypeChargeMap = grapheneReaderXYZPDB.getChargeCoeff();
        coeffPotentialMapNew = grapheneReaderXYZPDB.getCoeffPotential();
        bondTypeMap = grapheneReaderXYZPDB.getBondTypeMap();
        angleTypeMap = grapheneReaderXYZPDB.getAngleTypeMap();
        torsionTypeMap = grapheneReaderXYZPDB.getTorsionTypeMap();
        improperTypeMap = grapheneReaderXYZPDB.getImproperTypeMap();
        bondPotential = grapheneReaderXYZPDB.getBondPotential();
        anglePotential = grapheneReaderXYZPDB.getAnglePotential();
        torsionPotential = grapheneReaderXYZPDB.getTorsionPotential();
        improperPotential = grapheneReaderXYZPDB.getImproperPotential();
        uniqueElements = grapheneReaderXYZPDB.getUniqueElements();
        uniqueAtomTypeMap = grapheneReaderXYZPDB.getUniqueAtomTypeMap();

        newAtomTypeMap = grapheneReaderXYZPDB.getnewAtomTypeMap();
        typeChargeMap = grapheneReaderXYZPDB.getTypeChargeMap();
        Map<Integer, Double> chargeMap = grapheneReaderXYZPDB.getChargeMap();
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
         //   System.out.println(position +  " " + typeNew + "  " + i);
            speciesBuilderNew.addAtom(typeNew, position,  "");
        }
        species= speciesBuilderNew.setDynamic(true).build();
        setBondTypeMap(bondTypeMap);
        setAngleTypeMap(angleTypeMap);
        setTorsionTypeMap(torsionTypeMap);
        setImproperTypeMap(improperTypeMap);
        setUniqueAtomTypeMap(uniqueAtomTypeMap);
        setPotentials(bondPotential, anglePotential,torsionPotential, improperPotential,elementNumMap,positionMap );
        setNewMaps(atomMapGraphene, coeffPotentialMapNew, newAtomTypeMap, newTypeChargeMap, typeChargeMap);
        return species;
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


    public void makeBondingPotential(GrapheneReaderXYZPDB grapheneReaderXYZPDB, ISpecies species1, PotentialMasterBonding pmBonding, boolean ifAmberData){

        Map<Integer, List<int[]>> bondTypesMap2 = grapheneReaderXYZPDB.getBondTypeMap();
        Map<Integer, List<int[]>> angleTypesMap2 = grapheneReaderXYZPDB.getAngleTypeMap();
        Map<Integer, List<int[]>> torsionTypesMap2 = grapheneReaderXYZPDB.getTorsionTypeMap();
        Map<Integer, List<int[]>> improperTypesMap2 = grapheneReaderXYZPDB.getImproperTypeMap();
        Map<Integer, double[]> bondPotential = grapheneReaderXYZPDB.getBondPotential();
        Map<Integer, double[]> anglePotential = grapheneReaderXYZPDB.getAnglePotential();
        Map<Integer, double[]> torsionPotential = grapheneReaderXYZPDB.getTorsionPotential();
        Map<Integer, double[]> improperPotential = grapheneReaderXYZPDB.getImproperPotential();
        Unit kcals_mol = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        for (int i = 1; i < bondTypesMap2.size()+1; i++){
            double[] potential = bondPotential.get(i);

            P2HarmonicUFF p2HarmonicUFF = new P2HarmonicUFF(kcals_mol.toSim(potential[0]), Angstrom.UNIT.toSim(potential[1]));
            pmBonding.setBondingPotentialPair(species1, p2HarmonicUFF, bondTypesMap2.get(i));
        }

      /*  for (int i = 1; i < angleTypesMap2.size()+1; i++){
            double[] potential = anglePotential.get(i);
            P3BondAngle p3BondAngle = new P3BondAngle(Degree.UNIT.toSim(potential[1]), kcals_mol.toSim(potential[0]) );
            pmBonding.setBondingPotentialTriplet(species1, p3BondAngle, angleTypesMap2.get(i));
        }
        for (int i =1; i < torsionTypesMap2.size()+1; i++){
            double[] potential = torsionPotential.get(i);
            if (potential.length == 4){
                P4BondTorsion1 p4BondTorsion1 = new P4BondTorsion1(kcals_mol.toSim(potential[1]), (int) potential[2], Degree.UNIT.toSim(potential[3]) );
                pmBonding.setBondingPotentialQuad(species1, p4BondTorsion1, torsionTypesMap2.get(i));
            } else if (potential.length == 7) {
                P4BondTorsion2 p4BondTorsion2 = new P4BondTorsion2(kcals_mol.toSim(potential[1]), (int) potential[2], Degree.UNIT.toSim(potential[3]), kcals_mol.toSim(potential[4]), (int) potential[5], Degree.UNIT.toSim(potential[6]) );
                pmBonding.setBondingPotentialQuad(species1, p4BondTorsion2, torsionTypesMap2.get(i));
            } else if (potential.length == 10) {
                P4BondTorsion3 p4BondTorsion3 = new P4BondTorsion3(kcals_mol.toSim(potential[1]), (int) potential[2], Degree.UNIT.toSim(potential[3]), kcals_mol.toSim(potential[4]), (int) potential[5], Degree.UNIT.toSim(potential[6]), kcals_mol.toSim(potential[7]), (int) potential[8], Degree.UNIT.toSim(potential[9]));
                pmBonding.setBondingPotentialQuad(species1, p4BondTorsion3, torsionTypesMap2.get(i));
            } else  {
               throw new RuntimeException("error torsion type");
            }
        }*/
    }


    public void makeBondingPotential(GrapheneReaderXYZPDB grapheneReaderXYZPDB, ISpecies species1, PotentialMasterBonding pmBonding){
        Map<Integer, List<int[]>> bondTypesMap2 = grapheneReaderXYZPDB.getBondTypeMap();
        Map<Integer, List<int[]>> angleTypesMap2 = grapheneReaderXYZPDB.getAngleTypeMap();
        Map<Integer, List<int[]>> torsionTypesMap2 = grapheneReaderXYZPDB.getTorsionTypeMap();
        Map<Integer, List<int[]>> improperTypesMap2 = grapheneReaderXYZPDB.getImproperTypeMap();
        Map<Integer, double[]> bondPotential = grapheneReaderXYZPDB.getBondPotential();
        Map<Integer, double[]> anglePotential = grapheneReaderXYZPDB.getAnglePotential();
        Map<Integer, double[]> torsionPotential = grapheneReaderXYZPDB.getTorsionPotential();
        Map<Integer, double[]> improperPotential = grapheneReaderXYZPDB.getImproperPotential();
        Unit kcals_mol = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        for (int i = 1; i < bondTypesMap2.size()+1; i++){
            double[] potential = bondPotential.get(i);

            P2HarmonicUFF p2HarmonicUFF = new P2HarmonicUFF(kcals_mol.toSim(potential[0]), potential[1]);
            pmBonding.setBondingPotentialPair(species1, p2HarmonicUFF, bondTypesMap2.get(i));
        }

        for (int i = 1; i < angleTypesMap2.size()+1; i++){
            double[] potential = anglePotential.get(i);
            P3BondAngle p3BondAngle = new P3BondAngle(Degree.UNIT.toSim(potential[1]), kcals_mol.toSim(potential[0]) );
            pmBonding.setBondingPotentialTriplet(species1, p3BondAngle, angleTypesMap2.get(i));
        }
        for (int i =1; i < torsionTypesMap2.size()+1; i++){
            double[] potential = torsionPotential.get(i);
            if (potential.length == 4){
                BondTorsionOPLS bondTorsionOPLS = new BondTorsionOPLS(kcals_mol.toSim(potential[0]), kcals_mol.toSim(potential[1]), kcals_mol.toSim(potential[2]), kcals_mol.toSim(potential[3]));
                pmBonding.setBondingPotentialQuad(species1, bondTorsionOPLS, torsionTypesMap2.get(i));
            } else  {
                throw new RuntimeException("error torsion type");
            }
        }
    }

    public void makeNBPotential(GrapheneReaderXYZPDB grapheneReaderXYZPDB, List<List<AtomType>> pairsAtoms, PotentialMaster potentialMaster){
        Map<Integer, String> uniqueAtomMap = grapheneReaderXYZPDB.getUniqueAtomTypeMap();
        Map<String, double[]> coeffPotential = grapheneReaderXYZPDB.getCoeffPotential();
        int uniqueAtomsNum = uniqueAtomMap.size();
        int combinations = uniqueAtomsNum * uniqueAtomsNum;
        double sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey;
        int k = 0;
        double truncatedRadius = 12.5;
        UFF uff = new UFF();
        LJUFF[] p2LJ = new LJUFF[combinations];
        P2Electrostatic[] P2Electrostatics = new P2Electrostatic[combinations];
        IPotential2[] p2lj = new IPotential2[combinations];
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT),Mole.UNIT);
       // P2SoftSphericalSumTruncated[] p2Trunc = new P2SoftSphericalSumTruncated[uniqueAtomsNum*uniqueAtomsNum];
        TruncationFactory tf = new TruncationFactoryForceShift(truncatedRadius);
        for(List<AtomType>individualPair: pairsAtoms){
            AtomType atomNameOne = individualPair.get(0);
            AtomType atomNameTwo = individualPair.get(1);
            String atomTypeStringOne = String.valueOf(atomNameOne);
            String atomTypeStringTwo = String.valueOf(atomNameTwo);
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);

            double[] iKey = coeffPotential.get(atomTypeOne);
            double[] jKey = coeffPotential.get(atomTypeTwo);
            epsilonIKey = kcals.toSim(iKey[1]);
            sigmaIKey = iKey[0];
            epsilonJKey = kcals.toSim(jKey[0]);
            sigmaJKey = jKey[1];

            p2LJ[k] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);


            p2lj[k] = tf.make(p2LJ[k]);

            potentialMaster.setPairPotential(atomNameOne, atomNameTwo, p2lj[k], new double[]{1, 0, 0, 1});
            k ++;
        }
    }

    public void makeNBElectroPotential(GrapheneReaderXYZPDB grapheneReaderXYZPDB, List<List<AtomType>> pairsAtoms, PotentialMaster potentialMaster, PotentialComputeEwaldFourier ewaldFourier, List<AtomType> atomTypes, PotentialComputeEwaldFourier.EwaldParams params){
        Map<Integer, String> uniqueAtomMap = grapheneReaderXYZPDB.getUniqueAtomTypeMap();
        Map<String, double[]> coeffPotential = grapheneReaderXYZPDB.getCoeffPotential();
        Map<String, Double> electroCharge = grapheneReaderXYZPDB.getChargeCoeff();
        int uniqueAtomsNum = uniqueAtomMap.size();
        int combinations = uniqueAtomsNum * uniqueAtomsNum;
        double sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey;
        int k = 0;
        double truncatedRadius = 12.5;
/*
        for (int i = 0; i < uniqueAtomMap.size(); i ++){
            String atmName = String.valueOf(atomTypes.get(i));
            int startIndex = atmName.indexOf("[") + 1;
            int endIndex = atmName.indexOf("]");
            String nameNew = atmName.substring(startIndex, endIndex);
            ewaldFourier.setCharge(atomTypes.get(i), Electron.UNIT.toSim(electroCharge.get(nameNew)));
        }
        ewaldFourier.setkCut(params.kCut);
        ewaldFourier.setAlpha(params.alpha);*/

        TruncationFactory tfLJ = new TruncationFactoryForceShift(truncatedRadius);

        UFF uff = new UFF();
        LJUFF[] p2LJ = new LJUFF[combinations];
        IPotential2[] p2Electrostatics = new IPotential2[combinations];
        IPotential2[] p2lj = new IPotential2[combinations];
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT),Mole.UNIT);

        for(List<AtomType>individualPair: pairsAtoms){
            AtomType atomNameOne = individualPair.get(0);
            AtomType atomNameTwo = individualPair.get(1);
            String atomTypeStringOne = String.valueOf(atomNameOne);
            String atomTypeStringTwo = String.valueOf(atomNameTwo);
            String atomTypeOne = atomTypeStringOne.substring(9, atomTypeStringOne.length() - 1);
            String atomTypeTwo = atomTypeStringTwo.substring(9, atomTypeStringTwo.length() - 1);

            double[] iKey = coeffPotential.get(atomTypeOne);
            double[] jKey = coeffPotential.get(atomTypeTwo);
            epsilonIKey = kcals.toSim(iKey[0]);
            sigmaIKey = iKey[1];
            epsilonJKey = kcals.toSim(jKey[0]);
            sigmaJKey = jKey[1];
            p2LJ[k] = uff.vdw(sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey);
          //  System.out.println(atomTypeOne + " "+ Arrays.toString(iKey));
         //   System.out.println(atomTypeTwo + " "+ Arrays.toString(jKey));
            double chargeA = electroCharge.get(atomTypeOne);
            double chargeB = electroCharge.get(atomTypeTwo);
          //  p2Electrostatics[k] = new P2Ewald1Real(Electron.UNIT.toSim(chargeA) * Electron.UNIT.toSim(chargeB), params.alpha);
          //  p2lj[k] = tfLJ.make(p2LJ[k], p2Electrostatics[k]);
            p2lj[k] = tfLJ.make(p2LJ[k]);
            potentialMaster.setPairPotential(atomNameOne, atomNameTwo, p2lj[k], new double[]{1, 0, 0, 1});
            k ++;
        }
    }


    public AtomType atomFactory(String atomType){
        AtomType atomType1 = new AtomType(Carbon.INSTANCE);
        return atomType1;
    }
    private void setBondTypeMap(Map<Integer, List<int[]>> bondTypeMap){this.bondTypeMap = bondTypeMap;}
    private void setAngleTypeMap(Map<Integer, List<int[]>> angleTypeMap){this.angleTypeMap = angleTypeMap;}
    private void setTorsionTypeMap(Map<Integer, List<int[]>> torsionTypeMap){this.torsionTypeMap = torsionTypeMap;}
    private void setImproperTypeMap(Map<Integer, List<int[]>> improperTypeMap){this.improperTypeMap = improperTypeMap;}
    public void setUniqueAtomTypeMap(Map<Integer, String> map) {this.uniqueAtomTypeMap = map;}
    private void setIntArr(IntArrayList[] carboxyIntArr, IntArrayList[] hydroxyIntArr, IntArrayList[] epoxyIntArr){
        this.epoxyIntArr = epoxyIntArr;
        this.hydroxyIntArr = hydroxyIntArr;
        this.carboxyIntArr = carboxyIntArr;
    }
    private void totalCharge(Map<Integer, Double> totalCharge){
        double avCharge = 0.0;
        for (double q : totalCharge.values()) {
            avCharge += q;
        }
        System.out.println("Total charge = " + avCharge);
    }

    public Map<String, Map<Double, String>> getTypeChargeMap(){return  typeChargeMap;}
    public Map<Integer, double[]> getBondPotential(){return bondPotential;}
    public Map<Integer, double[]> getAnglePotential(){return anglePotential;}
    public Map<Integer, double[]> getTorsionPotential(){return torsionPotential;}
    public Map<Integer, double[]> getImproperPotential(){return improperPotential;}
    public Map<Integer, List<int[]>> getBondTypeMap(){return bondTypeMap;}
    public Map<Integer, List<int[]>> getAngleTypeMap(){return angleTypeMap;}
    public Map<Integer, List<int[]>> getTorsionTypeMap(){return torsionTypeMap;}
    public Map<Integer, List<int[]>> getImproperTypeMap(){return improperTypeMap;}
    public Map<Integer, Double> getChargeMap(){return chargeMap;}
    public Map<Integer, String> getUniqueAtomTypeMap() {return uniqueAtomTypeMap;}
    public Map<String, Double> getChargeCoeff(){return newTypeChargeMap;}
    public Map<String, double[]> getCoeffPotential(){return coeffPotentialMapNew;}

    public Map<String, List<Integer>> getnewAtomTypeMap(){return newAtomTypeMap;}

}
