package etomica.potential.OPLS_AA;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.IElement;
import etomica.chem.elements.Oxygen;
import etomica.potential.UFF.PDBReaderReplica;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;

import java.util.*;

public class PDBDataExtracterOPLS {
    ISpecies species;
    public static Map<Integer, AtomType> atomIdentifierMap = new HashMap<>();
    Map<Integer, String> modifiedAtomIdentifierMap = new HashMap<>();
    public static ArrayList<Integer> bondList = new ArrayList<>();
    public ArrayList<Integer> bondsNum = new ArrayList<>();
    Map<Integer, String> modifiedAtomIdentifierMapNB = new HashMap<>();
    Map<Integer, Vector> positions = new HashMap<>();
    static Map<Integer, Double[]> vdwMap = new HashMap<>();
    Map<Integer, Double> chargeMap = new HashMap<>();
    public Map<Integer,String > atomMap= new HashMap<>();
    public Map<Integer, String> atomMapModified = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectivity = new ArrayList<>();
    public ArrayList<ArrayList<Integer>> connectivityModified = new ArrayList<>();

    public Map<Integer, double[]> getbondingPotMap (Map<Integer, ArrayList<Integer[]>> bondTypeMap, ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomIdentifierMap){
        Map<Integer, double[]> bondPotMap = new HashMap<>();
        for (int i=0; i<bondTypeMap.size(); i++){
            ArrayList<Integer[]> value = bondTypeMap.get(i);
            Integer[] bondZero = value.get(0);
            double[] bondsArray = new double[2];
            String atomOne = atomIdentifierMap.get(bondZero[0]);
            String atomTwo = atomIdentifierMap.get(bondZero[1]);
           /* if (atomOne.equals("C_3") && atomTwo.equals("H")){
                bondsArray = new double[]{340, 1.090};
            } else if (atomOne.equals("C_3") && atomTwo.equals("C_3")) {
                bondsArray = new double[]{268, 1.529};
            } else if (atomOne.equals("C_3") && atomTwo.equals("C_2")) {
                bondsArray = new double[]{549, 1.34};
            }else if (atomOne.equals("C_2") && atomTwo.equals("C_2")) {
                bondsArray = new double[]{385, 1.46};
            } else if (atomOne.equals("C_3") && atomTwo.equals("O_3")) {
                bondsArray = new double[]{320, 1.41};
            }else if (atomOne.equals("O_3") && atomTwo.equals("H")) {
                bondsArray = new double[]{553, 0.9450};
            }else if (atomOne.equals("C_Ar") && atomTwo.equals("H")) {
                bondsArray = new double[]{367,1.08};
            }else if (atomOne.equals("C_Ar") && atomTwo.equals("C_Ar")) {
                bondsArray = new double[]{469,1.40};
            } else if (atomOne.equals("H") && atomTwo.equals("H")) {
                bondsArray = new double[]{}
            }*/
            if (atomOne.equals("Ha") && atomTwo.equals("Ha")) {
                bondsArray = new double[]{2093, 1.10};
            } else if (atomOne.equals("C3") && atomTwo.equals("C3")) {
                bondsArray = new double[]{1269.019,1.535};
            } else if (atomOne.equals("C3") && atomTwo.equals("HC")) {
                bondsArray = new double[]{1412.20764,1.092};
            }else if (atomOne.equals("C3") && atomTwo.equals("OH")) {
                bondsArray = new double[]{1315.0738,1.426};
            }else if (atomOne.equals("C3") && atomTwo.equals("H1")) {
                bondsArray = new double[]{1406.3461,1.093};
            }else if (atomOne.equals("OH") && atomTwo.equals("HO")) {
                bondsArray = new double[]{1547.4412,0.9740};
            }else if (atomOne.equals("N2") && atomTwo.equals("O")) {
                bondsArray = new double[]{3307.1533, 1.209};
            }else if (atomOne.equals("N2") && atomTwo.equals("OH")) {
                bondsArray = new double[]{1742.54616, 1.394};
            }else if (atomOne.equals("N1") && atomTwo.equals("N1")) {
                bondsArray = new double[]{5115.01356,1.124};
            }else if (atomOne.equals("C") && atomTwo.equals("O")) {
                bondsArray = new double[]{2713.0464,1.214};
            }else if (atomOne.equals("O") && atomTwo.equals("O")) {
                bondsArray = new double[]{1608.9872,1.43};
            }else if (atomOne.equals("C1") && atomTwo.equals("O")) {
                bondsArray = new double[]{3253.1436,1.166};
            }else if (atomOne.equals("N3") && atomTwo.equals("HN")) {
                bondsArray = new double[]{1650.1788,1.018};
            }else if (atomOne.equals("OW") && atomTwo.equals("HO")) {
                bondsArray = new double[]{2093.40,1.10};
            }else if (atomOne.equals("CA") && atomTwo.equals("CA")) {
                bondsArray = new double[]{2002.96512,1.387};
            }else if (atomOne.equals("CA") && atomTwo.equals("HA")) {
                bondsArray = new double[]{1441.51524,1.087};
            }else if (atomOne.equals("C2") && atomTwo.equals("C2")) {
                bondsArray = new double[]{2468.95596,1.324};
            }else if (atomOne.equals("C2") && atomTwo.equals("HA")) {
                bondsArray = new double[]{1441.51524,1.087};
            }else if (atomOne.equals("C3") && atomTwo.equals("C2")) {
                bondsArray = new double[]{1374.52644,1.508};
            }else if (atomOne.equals("C3") && atomTwo.equals("CA")) {
                bondsArray = new double[]{1354.4298,1.513};
            }else if (atomOne.equals("C2") && atomTwo.equals("CE")) {
                bondsArray = new double[]{2346.7014,1.339};
            }else if (atomOne.equals("CE") && atomTwo.equals("CE")) {
                bondsArray = new double[]{1634.9454,1.451};
            }else if (atomOne.equals("CE") && atomTwo.equals("HA")) {
                bondsArray = new double[]{1441.5152,1.087};
            }else if (atomOne.equals("") && atomTwo.equals("")) {
                bondsArray = new double[]{};
            }
            bondPotMap.put(i, bondsArray);
        }
        return bondPotMap;
    }

    public Map<Integer, double[]> getAngleMap (Map<Integer, ArrayList<Integer[]>> angleTypeMap, ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomIdentifierMap){
        Map<Integer, double[]> bondPotMap = new HashMap<>();
        double[] angleArray = new double[2];
        for (int i=0; i<angleTypeMap.size(); i++){
            ArrayList<Integer[]> value = angleTypeMap.get(i);
            Integer[] bondZero = value.get(0);
            String atomOne = atomIdentifierMap.get(bondZero[0]);
            String atomTwo = atomIdentifierMap.get(bondZero[1]);
            String atomThree = atomIdentifierMap.get(bondZero[2]);
            /*if(atomOne.equals("H") && atomTwo.equals("C_3") && atomThree.equals("H")){
                angleArray = new double[]{33,107.8};
            } else if (atomOne.equals("H") && atomTwo.equals("C_3") && atomThree.equals("C_3")) {
                angleArray = new double[]{37.8,110.7};
            } else if (atomOne.equals("C_3") && atomTwo.equals("C_3") && atomThree.equals("C_3")) {
                angleArray = new double[]{58.35, 112.7};
            } else if (atomOne.equals("H") && atomTwo.equals("C_2") && atomThree.equals("H")) {
                angleArray = new double[]{35,117};
            } else if (atomOne.equals("C_2") && atomTwo.equals("C_2") && atomThree.equals("H")) {
                angleArray = new double[]{35,120};
            }else if (atomOne.equals("C_3") && atomTwo.equals("C_2") && atomThree.equals("C_2")) {
                angleArray = new double[]{70,124};
            }else if (atomOne.equals("C_3") && atomTwo.equals("O_3") && atomThree.equals("H")) {
                angleArray = new double[]{55, 108};
            }else if (atomOne.equals("O_3") && atomTwo.equals("C_3") && atomThree.equals("H")) {
                angleArray = new double[]{35,109.5};
            }else if (atomOne.equals("H") && atomTwo.equals("C_Ar") && atomThree.equals("C_Ar")) {
                angleArray = new double[]{35,120};
            }else if (atomOne.equals("C_Ar") && atomTwo.equals("C_Ar") && atomThree.equals("C_Ar")) {
                angleArray = new double[]{63, 120};
            }else if (atomOne.equals("C_3") && atomTwo.equals("C_Ar") && atomThree.equals("C_Ar")) {
                angleArray = new double[]{70,120};
            }*/
            if (atomOne.equals("C3") && atomTwo.equals("C3") && atomThree.equals("HC")) {
                angleArray = new double[]{110.05,194.141916};
            }else if (atomOne.equals("HC") && atomTwo.equals("C3") && atomThree.equals("HC")) {
                angleArray = new double[]{108.35,165.085524};
            }else if (atomOne.equals("C3") && atomTwo.equals("C3") && atomThree.equals("OH")) {
                angleArray = new double[]{109.43,283.530096};
            }else if (atomOne.equals("C3") && atomTwo.equals("C3") && atomThree.equals("C3")) {
                angleArray = new double[]{110.63,264.647628};
            }else if (atomOne.equals("C3") && atomTwo.equals("C3") && atomThree.equals("H1")) {
                angleArray = new double[]{110.07,194.100};
            }else if (atomOne.equals("OH") && atomTwo.equals("C3") && atomThree.equals("H1")) {
                angleArray = new double[]{109.88,213.401196};
            }else if (atomOne.equals("C3") && atomTwo.equals("OH") && atomThree.equals("HO")) {
                angleArray = new double[]{108.16,197.1564};
            }else if (atomOne.equals("O") && atomTwo.equals("N2") && atomThree.equals("OH")) {
                angleArray = new double[]{112.15,316.354};
            }else if (atomOne.equals("HN") && atomTwo.equals("N3") && atomThree.equals("HN")) {
                angleArray = new double[]{107.13,172.9148};
            }else if (atomOne.equals("O") && atomTwo.equals("S2") && atomThree.equals("O")) {
                angleArray = new double[]{116.17,178.106472};
            }else if (atomOne.equals("HO") && atomTwo.equals("OW") && atomThree.equals("HO")) {
                angleArray = new double[]{120,0.08376};
            }else if (atomOne.equals("H1") && atomTwo.equals("C3") && atomThree.equals("H1")) {
                angleArray = new double[]{109.55,164.038824};
            }else if (atomOne.equals("CA") && atomTwo.equals("CA") && atomThree.equals("CA")) {
                angleArray = new double[]{119.97,281.2692240};
            }else if (atomOne.equals("CA") && atomTwo.equals("CA") && atomThree.equals("HA")) {
                angleArray = new double[]{120.01,202.892328};
            }else if (atomOne.equals("C2") && atomTwo.equals("C2") && atomThree.equals("HA")) {
                angleArray = new double[]{120.94,209.5074720};
            }else if (atomOne.equals("HA") && atomTwo.equals("C2") && atomThree.equals("HA")) {
                angleArray = new double[]{117.65,159.182136};
            }else if (atomOne.equals("C3") && atomTwo.equals("C2") && atomThree.equals("C2")) {
                angleArray = new double[]{123.42,269.3368440};
            }else if (atomOne.equals("C3") && atomTwo.equals("C2") && atomThree.equals("HA")) {
                angleArray = new double[]{191.16928,117.3};
            }else if (atomOne.equals("C2") && atomTwo.equals("C3") && atomThree.equals("HC")) {
                angleArray = new double[]{196.9052,110.49};
            }else if (atomOne.equals("CE") && atomTwo.equals("C2") && atomThree.equals("HA")) {
                angleArray = new double[]{121.19,207.539676};
            }else if (atomOne.equals("C2") && atomTwo.equals("CE") && atomThree.equals("CE")) {
                angleArray = new double[]{123.08,275.198364};
            }else if (atomOne.equals("CE") && atomTwo.equals("CE") && atomThree.equals("HA")) {
                angleArray = new double[]{115.9,198.873};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")) {
                angleArray = new double[]{};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")) {
                angleArray = new double[]{};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")) {
                angleArray = new double[]{};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")) {
                angleArray = new double[]{};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")) {
                angleArray = new double[]{};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")) {
                angleArray = new double[]{};
            }
        }
        return bondPotMap;
    }

    public Map<Integer, double[]> getTorsionMap (Map<Integer, ArrayList<Integer[]>> torsionTypeMap, ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomIdentifierMap){
        Map<Integer, double[]> bondPotMap = new HashMap<>();
        double[] angleArray = new double[2];
        for (int i=0; i<torsionTypeMap.size(); i++){
            ArrayList<Integer[]> value = torsionTypeMap.get(i);
            Integer[] bondZero = value.get(0);
            String atomOne = atomIdentifierMap.get(bondZero[0]);
            String atomTwo = atomIdentifierMap.get(bondZero[1]);
            String atomThree = atomIdentifierMap.get(bondZero[2]);
            String atomFour = atomIdentifierMap.get(bondZero[3]);
            /*if(atomOne.equals("H") && atomTwo.equals("C_3") && atomThree.equals("H")){
                angleArray = new double[]{33,107.8};
            } else if (atomOne.equals("H") && atomTwo.equals("C_3") && atomThree.equals("C_3")) {
                angleArray = new double[]{37.8,110.7};
            } else if (atomOne.equals("C_3") && atomTwo.equals("C_3") && atomThree.equals("C_3")) {
                angleArray = new double[]{58.35, 112.7};
            } else if (atomOne.equals("H") && atomTwo.equals("C_2") && atomThree.equals("H")) {
                angleArray = new double[]{35,117};
            } else if (atomOne.equals("C_2") && atomTwo.equals("C_2") && atomThree.equals("H")) {
                angleArray = new double[]{35,120};
            }else if (atomOne.equals("C_3") && atomTwo.equals("C_2") && atomThree.equals("C_2")) {
                angleArray = new double[]{70,124};
            }else if (atomOne.equals("C_3") && atomTwo.equals("O_3") && atomThree.equals("H")) {
                angleArray = new double[]{55, 108};
            }else if (atomOne.equals("O_3") && atomTwo.equals("C_3") && atomThree.equals("H")) {
                angleArray = new double[]{35,109.5};
            }else if (atomOne.equals("H") && atomTwo.equals("C_Ar") && atomThree.equals("C_Ar")) {
                angleArray = new double[]{35,120};
            }else if (atomOne.equals("C_Ar") && atomTwo.equals("C_Ar") && atomThree.equals("C_Ar")) {
                angleArray = new double[]{63, 120};
            }else if (atomOne.equals("C_3") && atomTwo.equals("C_Ar") && atomThree.equals("C_Ar")) {
                angleArray = new double[]{70,120};
            }*/
            if (atomOne.equals("OH") && atomTwo.equals("C3") && atomThree.equals("C3") && atomFour.equals("HC")) {
                angleArray = new double[]{0};
            }else if (atomOne.equals("C3") && atomTwo.equals("C3") && atomThree.equals("C3")&& atomFour.equals("HC")) {
                angleArray = new double[]{0.66988};
            }else if (atomOne.equals("HC") && atomTwo.equals("C3") && atomThree.equals("C3")&& atomFour.equals("H1")) {
                angleArray = new double[]{0.65128};
            }else if (atomOne.equals("HC") && atomTwo.equals("C3") && atomThree.equals("C3")&& atomFour.equals("OH")) {
                angleArray = new double[]{0.000};
            }else if (atomOne.equals("C3") && atomTwo.equals("C3") && atomThree.equals("OH")&& atomFour.equals("HO")) {
                angleArray = new double[]{0.66988};
            }else if (atomOne.equals("H1") && atomTwo.equals("C3") && atomThree.equals("OH")&& atomFour.equals("HO")) {
                angleArray = new double[]{0.6978};
            }else if (atomOne.equals("HC") && atomTwo.equals("C3") && atomThree.equals("C3")&& atomFour.equals("HC")) {
                angleArray = new double[]{0.6280200};
            }else if (atomOne.equals("C3") && atomTwo.equals("C3") && atomThree.equals("C3")&& atomFour.equals("C3")) {
                angleArray = new double[]{0.75362400};
            }else if (atomOne.equals("CA") && atomTwo.equals("CA") && atomThree.equals("CA")&& atomFour.equals("CA")) {
                angleArray = new double[]{15.17715};
            }else if (atomOne.equals("CA") && atomTwo.equals("CA") && atomThree.equals("CA")&& atomFour.equals("HA")) {
                angleArray = new double[]{15.17715};
            }else if (atomOne.equals("HA") && atomTwo.equals("CA") && atomThree.equals("CA")&& atomFour.equals("HA")) {
                angleArray = new double[]{15.17715};
            }else if (atomOne.equals("HA") && atomTwo.equals("C2") && atomThree.equals("C2")&& atomFour.equals("HA")) {
                angleArray = new double[]{27.84222};
            }else if (atomOne.equals("C3") && atomTwo.equals("C2") && atomThree.equals("C2")&& atomFour.equals("HA")) {
                angleArray = new double[]{27.84222};
            }else if (atomOne.equals("C2") && atomTwo.equals("C2") && atomThree.equals("C3")&& atomFour.equals("HC")) {
                angleArray = new double[]{1.590984};
            }else if (atomOne.equals("HA") && atomTwo.equals("C2") && atomThree.equals("C3")&& atomFour.equals("HC")) {
                angleArray = new double[]{0.0000};
            }else if (atomOne.equals("C3") && atomTwo.equals("CA") && atomThree.equals("CA")&& atomFour.equals("HA")) {
                angleArray = new double[]{15.17715};
            }else if (atomOne.equals("CE") && atomTwo.equals("CE") && atomThree.equals("C2")&& atomFour.equals("HA")) {
                angleArray = new double[]{27.84222};
            }else if (atomOne.equals("HA") && atomTwo.equals("C2") && atomThree.equals("CE")&& atomFour.equals("HA")) {
                angleArray = new double[]{27.844};
            }else if (atomOne.equals("C2") && atomTwo.equals("CE") && atomThree.equals("CE")&& atomFour.equals("HA")) {
                angleArray = new double[]{4.186};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")&& atomFour.equals("")) {
                angleArray = new double[]{};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")&& atomFour.equals("")) {
                angleArray = new double[]{};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")&& atomFour.equals("")) {
                angleArray = new double[]{};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")&& atomFour.equals("")) {
                angleArray = new double[]{};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")&& atomFour.equals("")) {
                angleArray = new double[]{};
            }else if (atomOne.equals("") && atomTwo.equals("") && atomThree.equals("")&& atomFour.equals("")) {
                angleArray = new double[]{};
            }
        }
        return bondPotMap;
    }
    public Map<Integer, double[]> getLJPotMap (List<String> listUniqueTypes, ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomIdentifierMap){
        Map<Integer, double[]> bondPotMap = new HashMap<>();
        for (int i=0; i<listUniqueTypes.size(); i++){
            double[] ljArray = new double[2];
            switch (listUniqueTypes.get(i)){
                case "C_3":
                    ljArray = new double[]{0.066, 3.5};
                case "H":
                    ljArray = new double[]{0.030,2.5};
                case "C_2":
                    ljArray = new double[]{0.076, 3.55};
                case "H_O":
                    ljArray = new double[]{0,0};
                case "O_3":
                    ljArray = new double[]{0.17,3.12};
                case "H_Ar":
                    ljArray = new double[]{0.03,2.42};
                case "C_Ar":
                    ljArray = new double[]{0.07, 3.55};
            }
            bondPotMap.put(i, ljArray);
        }
        return bondPotMap;
    }

    public Map<Integer, Double> chargePropane (String confName){
        Map<Integer, Double> chargePropane = new HashMap<>();
        switch (confName){
            case "propane":
                chargePropane.put(0, -0.065526);
                chargePropane.put(1, -0.058766);
                chargePropane.put(2, 0.022969);
                chargePropane.put(3, 0.022969);
                chargePropane.put(4, 0.022969);
                chargePropane.put(5, 0.022969);
                chargePropane.put(6, 0.022969);
                chargePropane.put(7, 0.022969);
                chargePropane.put(8, -0.065526);
                chargePropane.put(9, 0.026001);
                chargePropane.put(10, 0.026001);
            case "no":
                chargePropane.put(0,0.0);
                chargePropane.put(1,0.0);
            case "no2":
                chargePropane.put(0,0.0);
                chargePropane.put(1,0.0);
                chargePropane.put(2,0.0);
            case "ch4":
                chargePropane.put(0,0.0);
                chargePropane.put(1,0.0);
                chargePropane.put(2,0.0);
                chargePropane.put(3,0.0);
            case "n2":
                chargePropane.put(0,0.0);
                chargePropane.put(1,0.0);
            case "so2":
                chargePropane.put(0,0.0);
                chargePropane.put(1,0.0);
                chargePropane.put(2,0.0);
            case "co":
                chargePropane.put(0,0.0);
                chargePropane.put(1,0.0);
            case "co2":
                chargePropane.put(0,0.0);
                chargePropane.put(1,0.0);
                chargePropane.put(2,0.0);
        }

        return chargePropane;
    }




    public Map<Integer, String> atomIdentifierMapModifiedMaker (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified, String confName){
        Map<Integer, AtomType> atomIdentifierMap = atomIdentifier(connectivityModified, atomMapModified, confName);
        System.out.println(atomIdentifierMap + " in pdbExtractor");
        System.exit(1);
        for (Map.Entry<Integer, AtomType> entry : atomIdentifierMap.entrySet()) {
            String value = entry.getValue().toString();
            value = value.replace("AtomType[", "").replace("]", "");
            modifiedAtomIdentifierMap.put(entry.getKey(), value);
        }
        System.out.println(modifiedAtomIdentifierMap +" in extractor");
        return modifiedAtomIdentifierMap;
    }
    public Map<Integer, String> atomIdentifierMapModifiedMakerforLJ(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        PDBDataExtracterOPLS PDBDataExtracterOPLS = new PDBDataExtracterOPLS();
        PDBDataExtracterOPLS.elementInsideForNonBonded(connectivityModified, atomMapModified);
        for (Map.Entry<Integer, AtomType> entry : atomIdentifierMap.entrySet()) {
            String value = entry.getValue().toString();
            value = value.replace("AtomType[", "").replace("]", "");
            modifiedAtomIdentifierMap.put(entry.getKey(), value);
        }
        return modifiedAtomIdentifierMap;
    }
    protected Map<Integer, AtomType> atomIdentifier( ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified, String confName){
        int counter =0;
        ArrayList <Integer> atomNumbers = new ArrayList<>();
        if(connectivityModified.size() == 0){
            String element = atomMapModified.get(1);
            if(element.equals("Ar")){
                IElement Argon = null;
                AtomType Ar = new AtomType(Argon, "Ar");
                //System.out.println("The atom " + (i) +" is C_3" );
                atomIdentifierMap.put(1, Ar);
            } else if (element.equals("Ne")) {
                IElement Neon = null;
                AtomType Ne = new AtomType(Neon, "Ne");
                //System.out.println("The atom " + (i) +" is C_3" );
                atomIdentifierMap.put(1, Ne);
            } else {
                IElement Helium = null;
                AtomType He = new AtomType(Helium, "He");
                //System.out.println("The atom " + (i) +" is C_3" );
                atomIdentifierMap.put(1, He);
            }
        }
        System.out.println(confName);
        String[] parts = confName.split("//");
        // Return the last element
        String nameMolecule = parts[3];
        System.out.println(parts[3]);
        boolean ifAromatic = false;
        if(confName.equals("ethylbenzene") || confName.equals("oxylene") || confName.equals("mxylene") || confName.equals("pxylene") || confName.equals("toluene") || confName.equals("benzene")){
            ifAromatic = true;
        }


        System.out.println(connectivityModified + " connectivity Modified in atomIdentifier");
        System.out.println(atomMapModified);
        for(int i=0; i<connectivityModified.size(); i++){
            Integer retriveArrayFirstElementNumber = connectivityModified.get(i).get(0);
            String retriveArrayFirstElementName = atomMapModified.get(retriveArrayFirstElementNumber);
            System.out.println(retriveArrayFirstElementName + " :retriveFirstName");
            elementInsideSingle(connectivityModified, atomMapModified, i, retriveArrayFirstElementName, ifAromatic);
        }

        return atomIdentifierMap;
    }
    public void elementInsideSingle( ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified, int i,String retriveArrayFirstElementName, boolean ifAromatic){
        ArrayList<Integer>arrayInside= connectivityModified.get(i);
        int numCarb = 0, numHydro =0, numNitro=0, numOxy=0;
       /* for(int j =1; j<arrayInside.size(); j++){
            int arrayInsideElement = arrayInside.get(j);
            String arrayInsideElementName = atomMapModified.get(arrayInsideElement);
            if(arrayInsideElementName.equals("C")){
                numCarb++;
            } else if (arrayInsideElementName.equals("O")) {
                numOxy++;
            }else if (arrayInsideElementName.equals("N")) {
                numNitro++;
            }else {
                numHydro ++;
            }
        }*/
        int arraySize = connectivityModified.get(i).size();
        if(retriveArrayFirstElementName.equals("C")){
            if(arraySize == 5){
                AtomType C3 = new  AtomType(Carbon.INSTANCE, "C3");
                //System.out.println("The atom " + (i)+" is C_1 "  );
                atomIdentifierMap.put(i, C3);
            }else if(arraySize == 4 ){
                if(ifAromatic){
                    AtomType CA = new  AtomType(Carbon.INSTANCE, "CA");
                    //System.out.println("The atom " + (i)+" is C_1 "  );
                    atomIdentifierMap.put(i, CA);

                } else{
                    numHydro=0;
                    int numOne = connectivityModified.get(i).get(1);
                    int numTwo = connectivityModified.get(i).get(2);
                    int numThree = connectivityModified.get(i).get(3);
                    String atomOne = atomMapModified.get(numOne);
                    String atomTwo = atomMapModified.get(numTwo);
                    String atomThree = atomMapModified.get(numThree);
                    if(atomOne.equals("H")){
                        numHydro++;
                    }
                    if(atomTwo.equals("H")){
                        numHydro++;
                    }
                    if(atomThree.equals("H")){
                        numHydro++;
                    }
                    //  System.out.println(" numHydro " + numHydro);
                    //  System.out.println(" Here " +atomOne +" "+atomTwo + " " +atomThree);
                    //   System.out.println("Size "  +sizeOne + " " + sizeTwo + " " +sizeThree);

                    if (numHydro ==1 && atomMapModified.size() >9){
                        AtomType CE = new  AtomType(Carbon.INSTANCE, "CE");
                        //System.out.println("The atom " + (i)+" is C_1 "  );
                        atomIdentifierMap.put(i, CE);
                    } else {
                        // System.out.println(" Else Hydro " +numHydro);
                        //   System.out.println(atomMapModified.get(numOne) + " " +atomMapModified.get(numTwo) + " " +atomMapModified.get(numThree));
                        AtomType C2 = new  AtomType(Carbon.INSTANCE, "C2");
                        //System.out.println("The atom " + (i)+" is C_1 "  );
                        atomIdentifierMap.put(i, C2);
                    }
                }
            }
        } else if (retriveArrayFirstElementName.equals("H")) {
            int numConnect = connectivityModified.get(i).get(1);
            if(atomMapModified.get(numConnect).equals("C") && connectivityModified.get(i).size() == 4){
                AtomType HA = new  AtomType(Hydrogen.INSTANCE, "HA");
                //System.out.println("The atom " + (i)+" is C_1 "  );
                atomIdentifierMap.put(i, HA);
            } else if(atomMapModified.get(numConnect).equals("C") ){
                AtomType HC = new  AtomType(Hydrogen.INSTANCE, "HC");
                //System.out.println("The atom " + (i)+" is C_1 "  );
                atomIdentifierMap.put(i, HC);
            }else if(atomMapModified.get(numConnect).equals("O")){
                AtomType HO = new  AtomType(Hydrogen.INSTANCE, "HO");
                //System.out.println("The atom " + (i)+" is C_1 "  );
                atomIdentifierMap.put(i, HO);
            }
        }else if (retriveArrayFirstElementName.equals("O")) {
            int numConnect = connectivityModified.get(i).get(1);
            int numConnectTwo = connectivityModified.get(i).get(2);
            if(atomMapModified.get(numConnect).equals("H") && atomMapModified.get(numConnectTwo).equals("H")){
                AtomType Ow = new  AtomType(Oxygen.INSTANCE, "Ow");
                //System.out.println("The atom " + (i)+" is C_1 "  );
                atomIdentifierMap.put(i, Ow);
            }else  if(atomMapModified.get(numConnect).equals("H") ||atomMapModified.get(numConnectTwo).equals("H")){
                AtomType OH = new  AtomType(Oxygen.INSTANCE, "OH");
                //System.out.println("The atom " + (i)+" is C_1 "  );
                atomIdentifierMap.put(i, OH);
            }
        }
        if(ifAromatic){

        }else {

        }
    }

    public void elementInsideForNonBonded( ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        for(int i=0; i<connectivityModified.size(); i++){
            int numCarb =0, numHydr =0;
            ArrayList<Integer> insideElements = connectivityModified.get(i);
            int zerothElementNum = insideElements.get(0);
            String zerothElementName = atomMapModified.get(zerothElementNum);
            for(int j = 1; j<connectivityModified.get(i).size(); j++){
                int insideElementNum = connectivityModified.get(i).get(j);
                String insideElementName = atomMapModified.get(insideElementNum);
                if(insideElementName.equals("C")){
                    numCarb++;
                } else if (insideElementName.equals("H")){
                    numHydr++;
                }
            }
            String atomType = null;
            if(zerothElementName.equals("C")){
                if(connectivityModified.get(i).size() == 5){
                    if(numCarb == 0){
                        atomType = "C_0";
                        System.out.println("The atom " + (i) + " is C0");
                        modifiedAtomIdentifierMapNB.put(i,  atomType);
                    } else if (numCarb == 1) {
                        atomType = "C_1";
                        System.out.println("The atom " + (i) + " is C1");
                        modifiedAtomIdentifierMapNB.put(i,  atomType);
                    }else if (numCarb == 2) {
                        atomType = "C_2";
                        System.out.println("The atom " + (i) + " is C2");
                        modifiedAtomIdentifierMapNB.put(i,  atomType);
                    }else if (numCarb == 3) {
                        atomType = "C_2";
                        System.out.println("The atom " + (i) + " is C3");
                        modifiedAtomIdentifierMapNB.put(i,  atomType);
                    } else if (numCarb == 4) {
                        atomType = "C_4";
                        System.out.println("The atom " + (i) + " is C4");
                        modifiedAtomIdentifierMapNB.put(i,  atomType);
                    }
                }
                if(connectivityModified.get(i).size() == 4){
                    //C_200 - R2C=     C210 - RHC=    C220 - H2C=
                    if(numCarb == 0){
                        atomType = "C_220";
                        System.out.println("The atom " + (i) + " is C_220");
                        modifiedAtomIdentifierMapNB.put(i,  atomType);
                    } else if (numCarb == 1) {
                        atomType = "C_210";
                        System.out.println("The atom " + (i) + " is C210");
                        modifiedAtomIdentifierMapNB.put(i,  atomType);
                    }else if (numCarb == 2) {
                        atomType = "C_200";
                        System.out.println("The atom " + (i) + " is C200");
                        modifiedAtomIdentifierMapNB.put(i,  atomType);
                    }
                }
            } else if (zerothElementName.equals("H")){
                int firstElementNum = insideElements.get(1);
                String firstElementName = modifiedAtomIdentifierMapNB.get(firstElementNum);
                System.out.println((firstElementName == null) + " Start");
                if(firstElementName == null){
                    atomType = "H_2";
                    System.out.println("The atom " + (i) + " is H_2");
                    modifiedAtomIdentifierMapNB.put(i,  atomType);
                } else {
                    if(firstElementName.equals("C_0") || firstElementName.equals("C_1") || firstElementName.equals("C_2")|| firstElementName.equals("C_3")){
                        atomType = "H_2";
                        System.out.println("The atom " + (i) + " is H_2");
                        modifiedAtomIdentifierMapNB.put(i,  atomType);
                    } else if (firstElementName.equals("C_220") || firstElementName.equals("C_210")){
                        atomType = "H_2";
                        System.out.println("The atom " + (i) + " is H_2");
                        modifiedAtomIdentifierMapNB.put(i,  atomType);
                    }
                }

            }


        }
        System.out.println(modifiedAtomIdentifierMapNB + " NB");
    }

    public int priorityMapGenerator(String atomType){
        Map<String,Integer> priorityMap = new HashMap<>();
        priorityMap.put("C_0", 1);
        priorityMap.put("C_1", 2);
        priorityMap.put("C_2", 3);
        priorityMap.put("C_220", 4);
        priorityMap.put("C_210", 5);
        priorityMap.put("C_200", 6);
        priorityMap.put("H_2", 7);
        priorityMap.put("H_", 7);
        priorityMap.put("CT", 0);
        priorityMap.put("CTo", 1);
        priorityMap.put("CTb", 2);
        priorityMap.put("CTn", 3);
        priorityMap.put("CT0", 4);
        priorityMap.put("CT1", 5);
        priorityMap.put("CT2", 6);
        priorityMap.put("C1", 7);
        priorityMap.put("C2",8);
        priorityMap.put("C3", 9);
        priorityMap.put("N3", 10);
        priorityMap.put("N", 11);
        priorityMap.put("HC", 12);
        priorityMap.put("H2C", 13);
        priorityMap.put("H", 14);
        return priorityMap.get(atomType);
    }
    public List<int[]> bondSorter(List<int[]> duplets, Map<Integer, String> atomIdentifierMapModified){
        //List<int[]> priorityList = new ArrayList<>();
        System.out.println(atomIdentifierMapModified +" in Duplets");
        System.out.println(Arrays.deepToString(duplets.toArray()));
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
    public Map<Integer, String> getAtomIdentifierMapModified(){
        return modifiedAtomIdentifierMap;
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

    public ArrayList<Integer> setBondList(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        ArrayList<Integer> bondList = singleDoubleBondIdentifier(connectivityModified, atomMapModified);
        return bondList;
    }

    public ArrayList<Integer> singleDoubleBondIdentifier(ArrayList<ArrayList<Integer>> connectivity, Map<Integer, String> atomMap) {
        System.out.println("Already visited");
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
                    bondRequired = valenceHydrogen - (atomArraySize - 1);
                    // System.out.println("The hydrogen atom " + i + " is UNSATISFIED with " + bondRequired + " bond required");
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
                    //   System.out.println("The oxygen atom " + i + " is satisfied with one bond");
                    bondList.add(1);
                } else if (valenceOxygen == atomArraySize ) {
                    bondList.add(2);
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
                System.out.println(atomArraySize);
                if(atomArraySize == 3){
                    bondList.add(4);
                }
                if (valenceCarbon == atomArraySize - 1) {
                    // System.out.println("The carbon atom " + i + " is satisfied with one bond");
                    bondList.add(1);
                } else {
                    bondRequired = valenceCarbon - (atomArraySize - 1);
                    //   System.out.println("The carbon atom " + i + " is UNSATISFIED with " + bondRequired + " bond required");
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
                //  System.out.println(bondList.get(i) + " bondlist element " + i   );
                if (bondList.get(i) ==0){
                    System.out.println((i+1) + " " + bondList.get(i) + "i ");
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

    public double[] bondValueSender(String atomOne, String atomTwo){
        double[] bondValue = new double[2];
        if(atomOne.equals("C1") && atomTwo.equals("C1")){
            bondValue = new double[]{549.000, 1.340};
        } else if (atomOne.equals("C1") && atomTwo.equals("H2C")) {
            bondValue = new double[]{340.000, 1.080};
        } else if (atomOne.equals("CT1") && atomTwo.equals("HC")) {
            bondValue = new double[]{ 1412.20764, 1.092};
        }else if (atomOne.equals("CT1") && atomTwo.equals("CT1")) {
            bondValue = new double[]{ 1269.01908, 1.535};
        }else if (atomOne.equals("CT") && atomTwo.equals("HC")) {
            bondValue = new double[]{ 1412.20764, 1.092};
        } else if (atomOne.equals("C3") && atomTwo.equals("C3")) {
            bondValue = new double[]{1269.019,1.535};
        } else if (atomOne.equals("C2") && atomTwo.equals("C2")) {
            bondValue = new double[]{2468.956,1.324};
        }else if (atomOne.equals("C2") && atomTwo.equals("C3")) {
            bondValue = new double[]{1374.526,1.508};
        }else if (atomOne.equals("C2") && atomTwo.equals("Ce")) {
            bondValue = new double[]{2346.7,1.339};
        }else if (atomOne.equals("Ce") && atomTwo.equals("Ce")) {
            bondValue = new double[]{1634.94,1.451};
        }else if (atomOne.equals("Ca") && atomTwo.equals("Ca")) {
            bondValue = new double[]{2002.965,1.387};
        }else if (atomOne.equals("Ca") && atomTwo.equals("C3")) {
            bondValue = new double[]{1354.43,1.513};
        }else if (atomOne.equals("C3") && atomTwo.equals("Oh")) {
            bondValue = new double[]{1315.073,1.426};
        }else if (atomOne.equals("C1") && atomTwo.equals("O")) {
            bondValue = new double[]{3253.14,1.1666};
        }else if (atomOne.equals("C") && atomTwo.equals("O")) {
            bondValue = new double[]{2713,1.214};
        }else if (atomOne.equals("C3") && atomTwo.equals("Hc")) {
            bondValue = new double[]{1412.208,1.092};
        }else if (atomOne.equals("C3") && atomTwo.equals("H1")) {
            bondValue = new double[]{1406.346,1.093};
        }else if (atomOne.equals("C2") && atomTwo.equals("Ha")) {
            bondValue = new double[]{1441.515,1.087};
        }else if (atomOne.equals("Ce") && atomTwo.equals("Ha")) {
            bondValue = new double[]{1429.792,1.089};
        }else if (atomOne.equals("Ow") && atomTwo.equals("Ho")) {
            bondValue = new double[]{2093.4,1.1};
        }else if (atomOne.equals("O") && atomTwo.equals("O")) {
            bondValue = new double[]{1608.98,1.43};
        }else if (atomOne.equals("N1") && atomTwo.equals("N1")) {
            bondValue = new double[]{5115.014,1.124};
        }else if (atomOne.equals("Ha") && atomTwo.equals("Ha")) {
            bondValue = new double[]{2093.4,1.1};
        }else if (atomOne.equals("N3") && atomTwo.equals("Hn")) {
            bondValue = new double[]{1650.017,1.018};
        }else if (atomOne.equals("N2") && atomTwo.equals("O")) {
            bondValue = new double[]{3307.15,1.209};
        }else if (atomOne.equals("N2") && atomTwo.equals("Oh")) {
            bondValue = new double[]{1742.546,1.3194};
        }else if (atomOne.equals("S2") && atomTwo.equals("O")) {
            bondValue = new double[]{1396.716,1.599};
        }else if (atomOne.equals("Oh") && atomTwo.equals("Ho")) {
            bondValue = new double[]{1547.44,0.974};
        }
        return bondValue;
    }

    public double[] angleValueSender(String atomOne, String atomTwo, String atomThree){
        double[] angleValue = new double[2];
        if((atomOne.equals("C1") && atomTwo.equals("C1") && atomThree.equals("H2C") )|| (atomOne.equals("H2C") && atomTwo.equals("C1") && atomThree.equals("C1"))){
            angleValue = new double[]{549.000, 1.340};
        } else if (atomOne.equals("H2C") && atomTwo.equals("C1") && atomThree.equals("H2C")) {
            angleValue = new double[]{340.000, 1.080};
        } else if (atomOne.equals("HC") && atomTwo.equals("CT1") && atomThree.equals("HC")) {
            angleValue = new double[]{ 165.085524,  108.350};
        }else if (atomOne.equals("CT1") && atomTwo.equals("CT1") && atomThree.equals("HC")) {
            angleValue = new double[]{ 194.141916, 110.050};
        } else if (atomOne.equals("HC") && atomTwo.equals("CT") && atomThree.equals("HC")){
            angleValue = new double[]{ 165.085524,  108.350};
        } else if (atomOne.equals("C3") && atomTwo.equals("C3") && atomThree.equals("Hc")) {
            angleValue = new double[]{194.1419,  110.05};
        }else if (atomOne.equals("C3") && atomTwo.equals("C3") && atomThree.equals("H1")) {
            angleValue = new double[]{194.1,  110.07};
        }else if (atomOne.equals("C2") && atomTwo.equals("C2") && atomThree.equals("Ha")) {
            angleValue = new double[]{209.507,  120.94};
        }else if (atomOne.equals("C3") && atomTwo.equals("C2") && atomThree.equals("Ha")) {
            angleValue = new double[]{191.169,  117.3};
        }else if (atomOne.equals("C3") && atomTwo.equals("C2") && atomThree.equals("Hc")) {
            angleValue = new double[]{196.905,  110.49};
        }else if (atomOne.equals("Ce") && atomTwo.equals("C2") && atomThree.equals("Ha")) {
            angleValue = new double[]{207.497,  121.1};
        }else if (atomOne.equals("C3") && atomTwo.equals("C3") && atomThree.equals("Oh")) {
            angleValue = new double[]{283.53,  109.43};
        }else if (atomOne.equals("Ca") && atomTwo.equals("Ca") && atomThree.equals("Ha")) {
            angleValue = new double[]{202.89,  120.01};
        }else if (atomOne.equals("Ca") && atomTwo.equals("C3") && atomThree.equals("Hc")) {
            angleValue = new double[]{196.612,  110.15};
        }else if (atomOne.equals("Hc") && atomTwo.equals("C3") && atomThree.equals("Hc")) {
            angleValue = new double[]{165.0855,  108.35};
        }else if (atomOne.equals("Ha") && atomTwo.equals("C2") && atomThree.equals("Ha")) {
            angleValue = new double[]{159.182,  117.65};
        }else if (atomOne.equals("H1") && atomTwo.equals("C3") && atomThree.equals("H1")) {
            angleValue = new double[]{164.038,  109.55};
        }else if (atomOne.equals("C3") && atomTwo.equals("C3") && atomThree.equals("C3")) {
            angleValue = new double[]{264.647,  110.63};
        }else if (atomOne.equals("C3") && atomTwo.equals("C2") && atomThree.equals("C2")) {
            angleValue = new double[]{269.336,  123.42};
        }else if (atomOne.equals("C2") && atomTwo.equals("C2") && atomThree.equals("Ce")) {
            angleValue = new double[]{272.198,  123.08};
        }else if (atomOne.equals("Ca") && atomTwo.equals("Ca") && atomThree.equals("Ca")) {
            angleValue = new double[]{281.269,  119.97};
        }else if (atomOne.equals("Ca") && atomTwo.equals("Ca") && atomThree.equals("C3")) {
            angleValue = new double[]{267.28,  120.63};
        }else if (atomOne.equals("Ca") && atomTwo.equals("C3") && atomThree.equals("C3")) {
            angleValue = new double[]{264.8151,  112.09};
        }else if (atomOne.equals("Ho") && atomTwo.equals("Ow") && atomThree.equals("Ho")) {
            angleValue = new double[]{0.0837,  120.63};
        }else if (atomOne.equals("O") && atomTwo.equals("C1") && atomThree.equals("O")) {
            angleValue = new double[]{290.01,  180};
        }else if (atomOne.equals("Hn") && atomTwo.equals("N3") && atomThree.equals("Hn")) {
            angleValue = new double[]{172.914,  107.13};
        }else if (atomOne.equals("Oh") && atomTwo.equals("N2") && atomThree.equals("O")) {
            angleValue = new double[]{316.354,  112.15};
        }else if (atomOne.equals("O") && atomTwo.equals("S2") && atomThree.equals("O")) {
            angleValue = new double[]{178.106,  116.17};
        }else if (atomOne.equals("Oh") && atomTwo.equals("C3") && atomThree.equals("H1")) {
            angleValue = new double[]{213.4011,  109.88};
        }else if (atomOne.equals("C3") && atomTwo.equals("Oh") && atomThree.equals("Ho")) {
            angleValue = new double[]{197.1564,  108.16};
        }
        return angleValue;
    }

    public void listofNBInteractions(){
      /*  Map<String, Integer[]> tempMap = new HashMap<>();
        tempMap.put("Oxygen", new Integer[]{2, 5});
        tempMap.put("Hydrogen", new Integer[]{1, 3});
        tempMap.put("Helium", new Integer[]{4});

// chargeMap where Integer values map to Double
        Map<Integer, Double> chargeMap = new HashMap<>();
        chargeMap.put(1, 1.1);
        chargeMap.put(2, 2.2);
        chargeMap.put(3, 3.3);
        chargeMap.put(4, 4.4);
        chargeMap.put(5, 5.5);

// Create the new map to store the results
        Map<String, Double[]> newMap = new HashMap<>();

// Iterate through each entry in tempMap
        for (Map.Entry<String, Integer[]> entry : tempMap.entrySet()) {
            String key = entry.getKey();
            Integer[] integerArray = entry.getValue();

            // Create a new Double[] to store the corresponding charge values
            Double[] doubleArray = new Double[integerArray.length];

            // Replace Integer with Double from chargeMap
            for (int i = 0; i < integerArray.length; i++) {
                Integer intKey = integerArray[i];
                if (chargeMap.containsKey(intKey)) {
                    doubleArray[i] = chargeMap.get(intKey); // Map the Integer to its corresponding Double
                }
            }

            // Put the result in newMap
            newMap.put(key, doubleArray);
        }
        tempMap.forEach((key, value) -> System.out.println(key + "=" + Arrays.toString(value)));
// Print the newMap
       /* for (Map.Entry<String, Double[]> entry : newMap.entrySet()) {
            System.out.println(entry.getKey() + " -> " + Arrays.toString(entry.getValue()));
        }
        newMap.forEach((key, value) -> System.out.println(key + "=" + Arrays.toString(value)));*/
        // Example Map with String as key and Double[] as value
        Map<String, Double[]> map = new HashMap<>();
        map.put("String_1", new Double[]{1.0, 1.0, 1.0});
        map.put("String_2", new Double[]{2.0, 3.0, 1.0});
        map.put("String_3", new Double[]{7.0, 8.0, 7.0});

        // Create lists to store the updates (new keys and values)
        List<Map.Entry<String, Double[]>> newEntries = new ArrayList<>();

        // Check for unequal values and prepare new entries
        for (Map.Entry<String, Double[]> entry : map.entrySet()) {
            String key = entry.getKey();
            Double[] values = entry.getValue();

            // Check if all elements in the Double[] are the same
            boolean allEqual = true;
            Double firstValue = values[0];

            // Track unequal values
            List<Double> unequalValues = new ArrayList<>();

            // Identify unequal values
            for (Double value : values) {
                if (!value.equals(firstValue)) {
                    unequalValues.add(value);
                    allEqual = false;
                }
            }

            // If not all elements are the same, process the new key and update the map
            if (!allEqual) {
                // Create the original key's new value (only the repeated values)
                List<Double> repeatedValuesList = new ArrayList<>();
                repeatedValuesList.add(firstValue);
                for (Double value : values) {
                    if (value.equals(firstValue)) {
                        repeatedValuesList.add(value);
                    }
                }

                Double[] repeatedValues = repeatedValuesList.toArray(new Double[0]);
                newEntries.add(new AbstractMap.SimpleEntry<>(key, repeatedValues));

                // Now, add a new key for each of the unequal values
                int suffix = 1;
                for (Double unequalValue : unequalValues) {
                    String newKey = key + "_" + suffix++;
                    newEntries.add(new AbstractMap.SimpleEntry<>(newKey, new Double[]{unequalValue}));
                }

                // Debugging outputs
                System.out.println("Original key: " + key + " updated with values: ");
                for (Double v : repeatedValues) {
                    System.out.print(v + " ");
                }
                System.out.println();

                for (int i = 0; i < unequalValues.size(); i++) {
                    System.out.println("Created new key: " + key + "_" + (i+1) + " with value: ");
                    System.out.print(unequalValues.get(i) + "\n");
                }
            } else {
                // If all values are the same, just add the original key
                newEntries.add(entry);
            }
        }

        // Now, add the new entries to the map after the iteration is complete
        for (Map.Entry<String, Double[]> newEntry : newEntries) {
            map.put(newEntry.getKey(), newEntry.getValue());
        }

        // Print the updated map
        System.out.println("\nFinal Updated Map:");
        for (Map.Entry<String, Double[]> entry : map.entrySet()) {
            System.out.print(entry.getKey() + ": ");
            for (Double value : entry.getValue()) {
                System.out.print(value + " ");
            }
            System.out.println();
        }

    }

    public double[] torsionValueSender(String atomOne, String atomTwo, String atomThree, String atomFour){
        //kj/mol
        double[] torsionValue = new double[4];
        if(atomTwo.equals("CT1") && atomThree.equals("CT1")){
            torsionValue = new double[]{9, 1.40, 0.00, 3.0};
        }

        return torsionValue;
    }
    public double[] torsionValueSenderOPLS(String atomOne, String atomTwo, String atomThree, String atomFour){
        //kcal/mol
        double[] torsionValue = new double[3];
        if(atomOne.equals("HC") && atomTwo.equals("CT1") && atomThree.equals("CT1")&& atomFour.equals("HC")){
            torsionValue = new double[]{0,0,0.318};
        } else if (atomOne.equals("HC") && atomTwo.equals("CT1") && atomThree.equals("CT1")&& atomFour.equals("CT1")) {
            torsionValue = new double[]{0,0,0.366};
        }else if (atomOne.equals("CT1") && atomTwo.equals("CT1") && atomThree.equals("CT1")&& atomFour.equals("CT1")) {
            torsionValue = new double[]{1.740,-0.157,0.279};
        }
        return torsionValue;
    }

    public void vdwChargeMaker (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> modifiedAtomIdentifierMap){
        for(int i =0; i<modifiedAtomIdentifierMap.size(); i++){
            String atomName = modifiedAtomIdentifierMap.get(i);
            Double[] vdw = getVDWGAFF(atomName);
            vdwMap.put(i, vdw);
       }
    }

    public Map<String, Double> getChargeName(String confName, String atomname){
        int lastIndex = confName.lastIndexOf("//");

        // Extract the file name (substring from the position of last "//" to the end of the string)
        String fileName = confName.substring(lastIndex + 2); // +2 to skip the "//"

        Map<String, Double> map = new HashMap<>();
       /* if(fileName.equals("1propanol")){
            map =getChargeMapPropane();
        } else if (fileName.equals("ethanol")) {
            map = getChargeMapEthanol();
        } else*/ if (fileName.equals("ch4")) {
            map = getChargeMapMethane();
        }
        return  map;
    }

    public Double[] getVDWGAFF(String atomName){
     Double[] newOne = new Double[2];
        Map<String, Double[]> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put("C_3", new Double[]{1.9080, 0.1094});
        integerDoubleMap.put("H", new Double[]{1.4870, 0.0157});
        newOne = integerDoubleMap.get(atomName);
     return newOne;
    }

    public Map<Integer, Double> getChargeMapAmmonia(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,  -0.3436);
        integerDoubleMap.put(1,  0.1145);
        integerDoubleMap.put(2,  0.1145);
        integerDoubleMap.put(3,  0.1145);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapCO(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,  0.2174);
        integerDoubleMap.put(1, -0.2174);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapCO2(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,  0.3721);
        integerDoubleMap.put(1, -0.1860);
        integerDoubleMap.put(2, -0.1860);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapO2(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,  0.0);
        integerDoubleMap.put(1, 0.0);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapN2(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,  0.0);
        integerDoubleMap.put(1, 0.0);
        return integerDoubleMap;
    }

    public Map<Integer, Double> getChargeMapNO(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,  0.0);
        integerDoubleMap.put(1, 0.0);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapNO2(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,  -0.1018);
        integerDoubleMap.put(1, -0.1018);
        integerDoubleMap.put(2, 0.2037);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapSO2(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,  0.25188);
        integerDoubleMap.put(1, -0.1259);
        integerDoubleMap.put(2, -0.1259);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMap(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(1,  -0.068144);
        integerDoubleMap.put(2,  -0.068144);
        integerDoubleMap.put(3,  0.022715);
        integerDoubleMap.put(4,  0.022715);
        integerDoubleMap.put(5,  0.022715);
        integerDoubleMap.put(6,  0.022715);
        integerDoubleMap.put(7,  0.022715);
        integerDoubleMap.put(8,  0.022715);

        return integerDoubleMap;
    }

    public Map<String, Double> getChargeMapMethane(){
        //GAFF
        Map<String, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put("C_3p",  -0.077596);
        integerDoubleMap.put("H_p",  0.019399);
       /* integerDoubleMap.put(2,  0.019399);
        integerDoubleMap.put(3,  0.019399);
        integerDoubleMap.put(4,  0.019399);*/
        return integerDoubleMap;
    }


    public Map<Integer, Double> getChargeMapEthane(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,  -0.068144);
        integerDoubleMap.put(1,  -0.068144);
        integerDoubleMap.put(2,  0.022715);
        integerDoubleMap.put(3,  0.022715);
        integerDoubleMap.put(4,  0.022715);
        integerDoubleMap.put(5,  0.022715);
        integerDoubleMap.put(6,  0.022715);
        integerDoubleMap.put(7,  0.022715);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapPropane(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,   -0.058766);
        integerDoubleMap.put(1,   -0.065526);
        integerDoubleMap.put(2,   0.022969);
        integerDoubleMap.put(3,   0.022969);
        integerDoubleMap.put(4,   0.022969);
        integerDoubleMap.put(5,   0.022969);
        integerDoubleMap.put(6,   0.022969);
        integerDoubleMap.put(7,   0.022969);
        integerDoubleMap.put(8,   -0.058766);
        integerDoubleMap.put(9,   0.022969);
        integerDoubleMap.put(10,   0.022969);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapButane(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,   -0.065285);
        integerDoubleMap.put(1,   -0.056160);
        integerDoubleMap.put(2,   0.022977);
        integerDoubleMap.put(3,   0.022977);
        integerDoubleMap.put(4,   0.022977);
        integerDoubleMap.put(5,   0.022977);
        integerDoubleMap.put(6,   0.022977);
        integerDoubleMap.put(7,   0.022977);
        integerDoubleMap.put(8,   -0.056160);
        integerDoubleMap.put(9,   0.026256);
        integerDoubleMap.put(10,   0.026256);
        integerDoubleMap.put(11,   -0.065285);
        integerDoubleMap.put(12,   0.026256);
        integerDoubleMap.put(13,   0.026256);
        return integerDoubleMap;
    }

    public Map<Integer, Double> getChargeMapBenzene(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,   0.06175);
        integerDoubleMap.put(1,   -0.06175);
        integerDoubleMap.put(2,   -0.06175);
        integerDoubleMap.put(3,   0.06175);
        integerDoubleMap.put(4,   -0.06175);
        integerDoubleMap.put(5,   0.06175);
        integerDoubleMap.put(6,   -0.06175);
        integerDoubleMap.put(7,   0.06175);
        integerDoubleMap.put(8,   -0.06175);
        integerDoubleMap.put(9,   0.06175);
        integerDoubleMap.put(10,   -0.06175);
        integerDoubleMap.put(11,   0.06175);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapButadiene(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,   -0.0987);
        integerDoubleMap.put(1,   -0.0694);
        integerDoubleMap.put(2,   -0.0694);
        integerDoubleMap.put(3,   0.0534);
        integerDoubleMap.put(4,   0.06116);
        integerDoubleMap.put(5,   -0.0987);
        integerDoubleMap.put(6,   0.0534);
        integerDoubleMap.put(7,   0.06116);
        integerDoubleMap.put(8,   0.05348);
        integerDoubleMap.put(9,   0.05348);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapEthanol(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,  0.0413);
        integerDoubleMap.put(1,  -0.3952);
        integerDoubleMap.put(2,  0.2093);
        integerDoubleMap.put(3,  0.02516);
        integerDoubleMap.put(4,  0.02516);
        integerDoubleMap.put(5,  0.02516);
        integerDoubleMap.put(6,  -0.04182);
        integerDoubleMap.put(7,  0.05542);
        integerDoubleMap.put(8,  0.05542);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapEthene(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,  -0.1058);
        integerDoubleMap.put(1, -0.1058);
        integerDoubleMap.put(2,  0.05291);
        integerDoubleMap.put(3,  0.05291);
        integerDoubleMap.put(4,  0.05291);
        integerDoubleMap.put(5,  0.05291);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapEthylBenz(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,   0.02328);
        integerDoubleMap.put(1,   -0.0476);
        integerDoubleMap.put(2,   -0.05858);
        integerDoubleMap.put(3,   0.06202);
        integerDoubleMap.put(4,   -0.0615);
        integerDoubleMap.put(5,   0.0617);
        integerDoubleMap.put(6,   -0.0617);
        integerDoubleMap.put(7,   0.0617);
        integerDoubleMap.put(8,   -0.0615);
        integerDoubleMap.put(9,   0.0617);
        integerDoubleMap.put(10,   -0.05858);
        integerDoubleMap.put(11,   0.06202);
        integerDoubleMap.put(12,   -0.03052);
        integerDoubleMap.put(13,   0.02328);
        integerDoubleMap.put(14,   0.02328);
        integerDoubleMap.put(15,   -0.06126);
        integerDoubleMap.put(16,   0.031051);
        integerDoubleMap.put(17,   0.031051);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapIsoBut(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,   -0.3893);
        integerDoubleMap.put(1,   0.0257);
        integerDoubleMap.put(2,   0.0257);
        integerDoubleMap.put(3,   -0.0366);
        integerDoubleMap.put(4,   0.0257);
        integerDoubleMap.put(5,   0.2101);
        integerDoubleMap.put(6,   0.0577);
        integerDoubleMap.put(7,   -0.0366);
        integerDoubleMap.put(8,   0.0257);
        integerDoubleMap.put(9,   0.0257);
        integerDoubleMap.put(10,   0.0257);
        integerDoubleMap.put(11,   0.0257);
        integerDoubleMap.put(12,   -0.0366);
        integerDoubleMap.put(13,   0.0257);
        integerDoubleMap.put(14,   0.0257);
        return integerDoubleMap;
    }

    public Map<String, Double> getChargeMapMethanol(){
        //GAFF
        Map<String, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put("C_3p",   -0.3982);
        integerDoubleMap.put("O_3p",   0.2090);
        integerDoubleMap.put("H_p",   0.05208);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapMxylene(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,   0.02328);
        integerDoubleMap.put(1,   -0.0476);
        integerDoubleMap.put(2,   -0.05858);
        integerDoubleMap.put(3,   0.06202);
        integerDoubleMap.put(4,   -0.0615);
        integerDoubleMap.put(5,   0.0617);
        integerDoubleMap.put(6,   -0.0617);
        integerDoubleMap.put(7,   0.0617);
        integerDoubleMap.put(8,   -0.0615);
        integerDoubleMap.put(9,   0.0617);
        integerDoubleMap.put(10,   -0.05858);
        integerDoubleMap.put(11,   0.06202);
        integerDoubleMap.put(12,   -0.03052);
        integerDoubleMap.put(13,   0.02328);
        integerDoubleMap.put(14,   0.02328);
        integerDoubleMap.put(15,   -0.06126);
        integerDoubleMap.put(16,   0.031051);
        integerDoubleMap.put(17,   0.031051);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapOxylene(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,   0.02328);
        integerDoubleMap.put(1,   -0.0476);
        integerDoubleMap.put(2,   -0.05858);
        integerDoubleMap.put(3,   0.06202);
        integerDoubleMap.put(4,   -0.0615);
        integerDoubleMap.put(5,   0.0617);
        integerDoubleMap.put(6,   -0.0617);
        integerDoubleMap.put(7,   0.0617);
        integerDoubleMap.put(8,   -0.05858);
        integerDoubleMap.put(9,   0.0617);
        integerDoubleMap.put(10,   -0.04753);
        integerDoubleMap.put(11,   0.027762);
        integerDoubleMap.put(12,   -0.039460);
        integerDoubleMap.put(13,   0.02776);
        integerDoubleMap.put(14,   0.02776);
        integerDoubleMap.put(15,   -0.03946);
        integerDoubleMap.put(16,   0.027762);
        integerDoubleMap.put(17,   0.027762);
        return integerDoubleMap;
    }
    public Map<Integer, Double> getChargeMapPxylene(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,   0.02328);
        integerDoubleMap.put(1,   -0.0476);
        integerDoubleMap.put(2,   -0.05858);
        integerDoubleMap.put(3,   0.06202);
        integerDoubleMap.put(4,   -0.0615);
        integerDoubleMap.put(5,   0.0617);
        integerDoubleMap.put(6,   -0.0617);
        integerDoubleMap.put(7,   0.0617);
        integerDoubleMap.put(8,   -0.05858);
        integerDoubleMap.put(9,   0.0617);
        integerDoubleMap.put(10,   -0.04753);
        integerDoubleMap.put(11,   0.027762);
        integerDoubleMap.put(12,   -0.039460);
        integerDoubleMap.put(13,   0.02776);
        integerDoubleMap.put(14,   0.02776);
        integerDoubleMap.put(15,   -0.03946);
        integerDoubleMap.put(16,   0.027762);
        integerDoubleMap.put(17,   0.027762);
        return integerDoubleMap;
    }

    public Map<Integer, Double> getChargeMapPropene(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,   -0.09420);
        integerDoubleMap.put(1,   -0.10289);
        integerDoubleMap.put(2,   0.02716);
        integerDoubleMap.put(3,   0.02716);
        integerDoubleMap.put(4,  0.05316);
        integerDoubleMap.put(5,   0.05316);
        integerDoubleMap.put(6,   -0.04707);
        integerDoubleMap.put(7,   0.02716);
        integerDoubleMap.put(8,   0.05635);
        return integerDoubleMap;
    }

    public Map<Integer, Double> getChargeMapToluene(){
        //GAFF
        Map<Integer, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put(0,   0.02775);
        integerDoubleMap.put(1,   -0.05446);
        integerDoubleMap.put(2,   -0.058837);
        integerDoubleMap.put(3,   0.062020);
        integerDoubleMap.put(4,  -0.062020);
        integerDoubleMap.put(5,   0.06176);
        integerDoubleMap.put(6,   -0.06176);
        integerDoubleMap.put(7,   0.06176);
        integerDoubleMap.put(8,   -0.06176);
        integerDoubleMap.put(9,   0.0615);
        integerDoubleMap.put(10,   -0.0588);
        integerDoubleMap.put(11,   0.0277);
        integerDoubleMap.put(12,   -0.0397);
        integerDoubleMap.put(13,   -0.027753);
        integerDoubleMap.put(14,   0.06176);
        return integerDoubleMap;
    }



    public Double [] atomicPot (String atomtype){
        HashMap<String, Double[]> atomicConstant = new HashMap<>();
        atomicConstant.put("C_3", new Double[]{-0.240,3.5,0.066});
        atomicConstant.put("C_1", new Double[]{-0.180,3.5,0.066});
        atomicConstant.put("H", new Double[]{0.060,2.5,0.030});
        atomicConstant.put("C_210",new Double[]{-0.230,3.5,0.0760 });
      /*  atomicConstant.put("R2C", new double[]{-0.180,3.5,0.066});
        atomicConstant.put("R3C", new double[]{-0.060,3.5,0.066});
        atomicConstant.put("R4C", new double[]{0.000,3.5,0.066});
        atomicConstant.put("RH", new double[]{0.060,2.5,0.030});
        atomicConstant.put("ArC", new double[]{-0.115,3.550,0.070}); //Aromatic Carbon
        atomicConstant.put("ArH", new double[]{0.115,2.420,0.030}); //Aromatic H
        atomicConstant.put("ArRC", new double[]{-0.065,3.5,0.066});//Carbon attached to AromaticGroup
        atomicConstant.put("ArRC2", new double[]{-0.005,3.550,0.076}); //Ethyl carbon to Aromatic
        atomicConstant.put("R2C2", new double[]{-0.115,3.5,0.076}); //Secondary Carb Double bonded
        atomicConstant.put("RHC2", new double[]{-0.230,3.5,0.076}); //Primary Carb Double Bonded
        atomicConstant.put("H2C2", new double[]{0.000,3.5,0.066}); //Ethyl Carb Double Bonded
        atomicConstant.put("RC2H", new double[]{0.000,3.5,0.030}); //Ethyl Carb Double Bonded HYDROGEN
        atomicConstant.put("RO", new double[]{0.000,3.5,0.17}); //Oxygen from alcohol
        atomicConstant.put("ROH", new double[]{0.000,3.5,0}); //hydrogen attached to O
        atomicConstant.put("RCHOH", new double[]{0.000,3.5,0.03});//hydrogen attached to C in alcohol
        atomicConstant.put("RCH2OH", new double[]{0.000,3.5,0.066}); //Primary Alco C
        atomicConstant.put("R2CHOH", new double[]{0.000,3.5,0.066}); //Secondary Alco C
        atomicConstant.put("R3COH", new double[]{0.000,3.5,0.066}); //Tertiary Alco C
        atomicConstant.put("Ar(C)OH", new double[]{0.000,3.5,0.070}); //Phenol Carbon
        atomicConstant.put("ArC(O)H", new double[]{-0.585,3.070,0.170}); //Tertiary Alco C
        atomicConstant.put("ArCO(H)", new double[]{0.435,0,0}); //Tertiary Alco C*/
        Double [] sample = atomicConstant.get(atomtype);
        return sample;
    }
    public double getCharge(String atom){
        Map<String, Double> integerDoubleMap = new HashMap<>();
        integerDoubleMap.put("C_3",  -0.068144);
        integerDoubleMap.put("H",  0.022715);
        return integerDoubleMap.get(atom);
    }

    public ISpecies getSpeciesNew (String confName){
        PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
        PDBDataExtracterOPLS PDBDataExtracterOPLS = new PDBDataExtracterOPLS();
        species = pdbReaderReplica.getSpecies(confName, true, new Vector3D(5.0,5.0,5.0),false);
        System.out.println("In simulation OPLS");
        Map<Integer, String> atomMap = pdbReaderReplica.atomMap;
        HashMap<Integer, String> atomMapModified = pdbReaderReplica.atomMapModified;
        ArrayList<ArrayList<Integer>> connectivity = pdbReaderReplica.connectivity;
        ArrayList<ArrayList<Integer>> connectivityModified = pdbReaderReplica.connectivityModified;
        Map<Integer, Vector> positionMap = pdbReaderReplica.getPositions();
        setPositions(positionMap);
        System.out.println(connectivity + " here OPLS");
        System.out.println(connectivityModified + " here OPLS");
        System.out.println(atomMap + " here OPLS");
        System.out.println(atomMapModified + " here OPLS");
        Map<Integer, String> atomIdentifierMapModified = PDBDataExtracterOPLS.atomIdentifierMapModifiedMaker(connectivityModified, atomMapModified, confName);
        System.out.println(atomIdentifierMapModified + " in OPLS Simulation");
        PDBDataExtracterOPLS.elementInsideForNonBonded(connectivityModified, atomMapModified);
        List<int[]> duplets = pdbReaderReplica.getBondedAtomList(connectivityModified);
        System.out.println(Arrays.deepToString(duplets.toArray())+ ": listOfBonds");
        System.out.println(atomIdentifierMapModified);
        List<int[]> triplets = pdbReaderReplica.getAngleList(connectivityModified);
        System.out.println(Arrays.deepToString(triplets.toArray())+ ": listOfAngleModified");
        List<int[]> quadruplets = pdbReaderReplica.getTorsionList(connectivity);
        System.out.println(Arrays.deepToString(quadruplets.toArray())+ " listOfTorsionModified");
        ArrayList<Integer> bondList = pdbReaderReplica.getBondList(connectivity, atomMap);
        List<int[]> dupletsSorted= PDBDataExtracterOPLS.bondSorter(duplets, atomIdentifierMapModified);
        System.out.println(Arrays.deepToString(dupletsSorted.toArray()) + " duplets Sorted");
        List<int[]>tripletsSorted=PDBDataExtracterOPLS.angleSorter(triplets, atomIdentifierMapModified);
        System.out.println(Arrays.deepToString(tripletsSorted.toArray()) + "triplets");
        Map<Integer, String> modifiedAtomIdentifierMap = PDBDataExtracterOPLS.getAtomIdentifierMapModified();
        List<int[]>quadrupletsSorted=PDBDataExtracterOPLS.torsionSorter(quadruplets, atomIdentifierMapModified);
        Map<Integer, String>modifiedAtomIdentifierMapNB = PDBDataExtracterOPLS.getModifiedAtomIdentifierMapNB();
        System.out.println(modifiedAtomIdentifierMapNB + " NonBonded");
        return species;
    }
    public void setPositions(Map<Integer, Vector> positions){this.positions = positions;}
    public Map<Integer, Vector> getPositions (){return positions;}

    public Map<Integer, String> getModifiedAtomIdentifierMapNB(){
        return modifiedAtomIdentifierMapNB;
    }

    public void bondingMoleculeMethanol(){}


    public static void main(String[] args) {
        PDBDataExtracterOPLS PDBDataExtracterOPLS = new PDBDataExtracterOPLS();
        String confName = "F://Avagadro//molecule_GAFF//benzene";
        PDBDataExtracterOPLS.listofNBInteractions();
        //PDBDataExtracterOPLS.getSpeciesNew(confName);
        /*String confName = "F://ethane";
        species = PDBReader.getSpeciesAnotherFF(confName);
        Map<Integer, String> atomMap = PDBReader.atomMap;
        HashMap<Integer, String> atomMapModified = PDBReader.atomMapModified;
        ArrayList<ArrayList<Integer>> connectivity = PDBReader.connectivity;
        ArrayList<ArrayList<Integer>> connectivityModified = PDBReader.connectivityModified;
        System.out.println(connectivity + " here OPLS");
        System.out.println(connectivityModified + " here OPLS");
        System.out.println(atomMap + " here OPLS");
        System.out.println(atomMapModified + " here OPLS");
        Map<Integer, String> atomIdentifierMapModified = atomIdentifierMapModifiedMaker(connectivityModified, atomMapModified);
        elementInsideForNonBonded(connectivityModified, atomMapModified);
        List<int[]> duplets = PDBReader.getBondedAtomList(connectivityModified);
        System.out.println(Arrays.deepToString(duplets.toArray())+ ": listOfBonds");
        System.out.println(atomIdentifierMapModified);
        List<int[]> triplets = PDBReader.getAngleList(connectivityModified);
        System.out.println(Arrays.deepToString(triplets.toArray())+ ": listOfAngleModified");
        List<int[]> quadruplets = PDBReader.getTorsionList(connectivity);
        System.out.println(Arrays.deepToString(quadruplets.toArray())+ " listOfTorsionModified");
        ArrayList<Integer> bondList = PDBReader.getBondList(connectivity, atomMap);
        List<int[]> dupletsSorted= bondSorter(duplets, atomIdentifierMapModified);
        System.out.println(Arrays.deepToString(dupletsSorted.toArray()) + " duplets Sorted");
        List<int[]>tripletsSorted=angleSorter(triplets, atomIdentifierMapModified);
        System.out.println(Arrays.deepToString(tripletsSorted.toArray()) + "triplets");
        List<int[]>quadrupletsSorted=torsionSorter(quadruplets, atomIdentifierMapModified);
        System.out.println(modifiedAtomIdentifierMapNB + " NonBonded");
        vdwChargeMaker(connectivityModified, atomIdentifierMapModified);
        chargeMap = getChargeMap();*/
        /*System.out.println(PDBDataExtracterOPLS.chargeMap);
        for (Map.Entry<Integer, Double[]> entry : vdwMap.entrySet()) {
            Integer key = entry.getKey();
            Double[] value = entry.getValue();
            System.out.println(key + " : " + Arrays.toString(value));
        }*/
    }
    public Map<Integer, String> makeMethanolTypes (){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(2, "oh");
        atomTypes.put(1, "c3");//
        atomTypes.put(3, "ho");//
        atomTypes.put(4, "h1");//4,5,6
        atomTypes.put(5, "h1");
        atomTypes.put(6, "h1");
        return atomTypes;
    }
    public Map<Integer, String> makeEthanolTypes (){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(2, "oh");
        atomTypes.put(4, "hc");//4,5,6
        atomTypes.put(5, "hc");
        atomTypes.put(6, "hc");
        atomTypes.put(1, "c3");//1,7
        atomTypes.put(3, "ho");//
        atomTypes.put(8, "h1");//8,9
        atomTypes.put(9, "h1");
        atomTypes.put(7, "c31");
        return atomTypes;
    }

    public Map<Integer, String> makeOnePropanolTypes (){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "oh");
        atomTypes.put(2, "hc");//2,3,5,11,12
        atomTypes.put(4, "c3");//4,7,10
        atomTypes.put(6, "ho");//
        atomTypes.put(8, "h1");//8,9
        atomTypes.put(9, "h1");
        atomTypes.put(11, "hc1");//new hc 11,12
        atomTypes.put(12, "hc1");
        atomTypes.put(7, "c31");//7
        atomTypes.put(10, "c32");//10
        return atomTypes;
    }


    public Map<Integer, String> makeIsoButanolTypes (){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "oh");
        atomTypes.put(2, "hc");//2,3,5,9,10,11,12,14,15
        atomTypes.put(3, "hc");
        atomTypes.put(5, "hc");
        atomTypes.put(9, "hc");
        atomTypes.put(10, "hc");
        atomTypes.put(14, "hc");
        atomTypes.put(15, "hc");
        atomTypes.put(11, "hc");
        atomTypes.put(12, "hc");
        atomTypes.put(4, "c3");//4,7,8,13
        atomTypes.put(8, "c3");
        atomTypes.put(13, "c3");
        atomTypes.put(6, "ho");//6
        atomTypes.put(7, "c31");
        return atomTypes;
    }

    public Map<Integer, String> makeTwoPropanolTypes (){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "oh");
        atomTypes.put(2, "hc");//3,5,9,10,11
        atomTypes.put(3, "hc");
        atomTypes.put(5, "hc");
        atomTypes.put(9, "hc");
        atomTypes.put(10, "hc");
        atomTypes.put(11, "hc");
        atomTypes.put(4, "c3");//7,8
        atomTypes.put(8, "c3");
        atomTypes.put(6, "ho");
        atomTypes.put(12, "h1");
        atomTypes.put(7, "c31");//diff
        return atomTypes;
    }

    public Map<Integer, String> makeBenzeneTypes (){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "ha"); //1,4,6,8,10,12
        atomTypes.put(4, "ha");
        atomTypes.put(6, "ha");
        atomTypes.put(8, "ha");
        atomTypes.put(10, "ha");
        atomTypes.put(12, "ha");
        atomTypes.put(2, "ca"); //2,3,5,7,9,11
        atomTypes.put(3, "ca");
        atomTypes.put(5, "ca");
        atomTypes.put(7, "ca");
        atomTypes.put(9, "ca");
        atomTypes.put(11, "ca");
        return atomTypes;
    }

    public Map<Integer, String> makeTolueneTypes (){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(4, "ha"); //4,6,8,10,15
        atomTypes.put(2, "ca"); //2,3,5,7,9,11
        atomTypes.put(13, "c3"); //13
        atomTypes.put(1, "hc"); //1,12,14
        atomTypes.put(12, "hc");
        atomTypes.put(14, "hc");
        atomTypes.put(10, "ha");
        atomTypes.put(6, "ca1");
        atomTypes.put(8, "ca1");
        atomTypes.put(15, "ca1");
        atomTypes.put(5, "ca2");
        atomTypes.put(7, "ca2");
        atomTypes.put(9, "ca2");
        atomTypes.put(3, "ca3");
        atomTypes.put(11, "ca3");
        atomTypes.put(2, "ca4");
        //4,10
       //6, 8, 15
        //5,7,9
        //3,11
        //2
        return atomTypes;
    }

    public Map<Integer, String> makeEthylBenzeneTypes (){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(4, "ha"); //4.6.8.10.12
        atomTypes.put(6, "ha");
        atomTypes.put(8, "ha");
        atomTypes.put(10, "ha");
        atomTypes.put(12, "ha");
        atomTypes.put(2, "ca"); //2,3,5,7,9,11,
        atomTypes.put(13, "c3");//13,16
        atomTypes.put(16, "c31");
        atomTypes.put(1, "hc");//1.14.15.17.18
        atomTypes.put(14, "hc");
        atomTypes.put(15, "hc");
        atomTypes.put(17, "hc1");
        atomTypes.put(18, "hc1");
        atomTypes.put(3, "ca1");
        atomTypes.put(11, "ca1");
        atomTypes.put(5, "ca2");
        atomTypes.put(7, "ca2");
        atomTypes.put(9, "ca2");
        //2
        //3.11
        //5.7.9
        //1.14.15
        //17.18
        //13.
        //16
        return atomTypes;
    }

    public Map<Integer, String> makeOXyleneTypes (){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(4, "ha"); //4,6,8,10,
        atomTypes.put(6, "ha");
        atomTypes.put(8, "ha");
        atomTypes.put(10, "ha");
        atomTypes.put(2, "ca"); //2,3,5,7,9,11
        atomTypes.put(3, "ca1");
        atomTypes.put( 5, "ca2");
        atomTypes.put( 7, "ca2");
        atomTypes.put( 9, "ca1");
        atomTypes.put( 11, "ca");
        atomTypes.put(1 , "hc");//1,12,14,15,17,18
        atomTypes.put( 12, "hc");
        atomTypes.put( 14, "hc");
        atomTypes.put(15, "hc");
        atomTypes.put( 17, "hc");
        atomTypes.put( 18, "hc");
        atomTypes.put(13, "c3");//13,16
        atomTypes.put(16, "c3");
        //2,11
        //3,9
        //5,7
        return atomTypes;
    }

    public Map<Integer, String> makeMXyleneTypes (){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(4, "ha"); //4,6,10,15
        atomTypes.put(6, "ha");
        atomTypes.put(10, "ha");
        atomTypes.put(15, "ha");
        atomTypes.put(2, "ca"); //2,3,5,7,9,11
        atomTypes.put(9, "ca");
        atomTypes.put(3, "ca");
        atomTypes.put(7, "ca");
        atomTypes.put(5, "ca");
        atomTypes.put(11, "ca");
        atomTypes.put(1, "hc");//1,8,12,14,17,18
        atomTypes.put(8, "hc");
        atomTypes.put(12, "hc");
        atomTypes.put(14, "hc");
        atomTypes.put(17, "hc");
        atomTypes.put(18, "hc");
        atomTypes.put(13, "c3");//13,16
        atomTypes.put(16, "c3");
        //2.9
        //3.7
        //5
        //11
        return atomTypes;
    }

    public Map<Integer, String> makePXyleneTypes (){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(4, "ha"); //4,6,10,18
        atomTypes.put(6, "ha");
        atomTypes.put(10, "ha");
        atomTypes.put(18, "ha");
        atomTypes.put(3, "ca"); //2,3,5,7,9,11
        atomTypes.put(5, "ca");
        atomTypes.put(9, "ca");
        atomTypes.put(11, "ca");
        atomTypes.put(2, "ca1");
        atomTypes.put(7, "ca1");
        atomTypes.put(1, "hc");//1,8,12,14,16,17
        atomTypes.put(8, "hc");
        atomTypes.put( 12, "hc");
        atomTypes.put( 14, "hc");
        atomTypes.put( 16, "hc");
        atomTypes.put(17, "hc");
        atomTypes.put(13, "c3");//13,15
        atomTypes.put(15, "c3");
        //3.5.9.11
        //2,7
        return atomTypes;
    }

    public Map<Integer, String> makeMethaneTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1,  "c3"); //1
        atomTypes.put( 2,  "hc"); //2.3.4.5
        atomTypes.put( 3,  "hc");
        atomTypes.put( 4,  "hc");
        atomTypes.put( 5,  "hc");
        return atomTypes;
    }

    public Map<Integer, String> makeEthaneTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1,  "c3");//1.2
        atomTypes.put(2,  "c3");
        atomTypes.put(3, "hc"); //3.4.5.6.7.8
        atomTypes.put(4, "hc");
        atomTypes.put(5, "hc");
        atomTypes.put(6, "hc");
        atomTypes.put(7, "hc");
        atomTypes.put(8, "hc");
        return atomTypes;
    }

    public Map<Integer, String> makePropaneTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1,  "c3");//1.2.9
        atomTypes.put(9,  "c3");
        atomTypes.put(2,  "c31");
        atomTypes.put(3, "hc"); //3.4.5.6.7.8.10.11
        atomTypes.put(4, "hc");
        atomTypes.put( 5, "hc");
        atomTypes.put(6 , "hc");
        atomTypes.put( 7, "hc");
        atomTypes.put(8 , "hc");
        atomTypes.put( 10, "hc");
        atomTypes.put(11, "hc");
        //1.9
        //2
        //10.11
        return atomTypes;
    }

    public Map<Integer, String> makeButaneTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1,  "c3");//1,2,10,12
        atomTypes.put(12,  "c3");
        atomTypes.put(2,  "c31");
        atomTypes.put(9,  "c31");
        atomTypes.put(3,  "hc"); // 3, 4 5,6,7,8,10,11,13,14
        atomTypes.put(4,  "hc");
        atomTypes.put(5,  "hc");
        atomTypes.put(6,  "hc");
        atomTypes.put(7,  "hc");
        atomTypes.put(8,  "hc");
        atomTypes.put(10,  "hc");
        atomTypes.put(11,  "hc");
        atomTypes.put(13,  "hc");
        atomTypes.put(14,  "hc");

        //2,9
        //1,12
        //3,4,5,6,7,8
        //10,11,13,14
        return atomTypes;
    }

    public Map<Integer, String> makeEtheneTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1,  "c2");//1.2
        atomTypes.put(2,  "c2");
        atomTypes.put(3, "ha"); //3.4.5.6
        atomTypes.put(4, "ha");
        atomTypes.put(5, "ha");
        atomTypes.put(6, "ha");
        return atomTypes;
    }

    public Map<Integer, String> makePropeneTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1,  "c2");//1.2
        atomTypes.put(2,  "c21");
        atomTypes.put(7,  "c3"); // 7
        atomTypes.put(3,  "hc"); //3,4,8
        atomTypes.put(4,  "hc");
        atomTypes.put(8,  "hc");
        atomTypes.put(5, "ha"); //5.6.9
        atomTypes.put(6, "ha");
        atomTypes.put(9, "ha");
        //1
        //2
        return atomTypes;
    }

    public Map<Integer, String> makeButadieneTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1,  "c2");//1.6
        atomTypes.put(6,  "c2");
        atomTypes.put(2,  "ce"); // 2.3
        atomTypes.put(3,  "ce");
        atomTypes.put(4,  "ha"); //4.5.7.8.9.10
        atomTypes.put(5,  "ha");
        atomTypes.put(7,  "ha");
        atomTypes.put(8,  "ha");
        atomTypes.put(9,  "ha");
        atomTypes.put(10,  "ha");
        return atomTypes;
    }

    public Map<Integer, String> makeCOTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "c3");
        atomTypes.put(2, "oh");
        return atomTypes;
    }

    public Map<Integer, String> makeCOTwoTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "c1");
        atomTypes.put(2, "o");
        atomTypes.put(3, "o");
        return atomTypes;
    }

    public Map<Integer, String> makeNTwoTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "n1");
        atomTypes.put(2, "n1");
        return atomTypes;
    }

    public Map<Integer, String> makeAmmoniaTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "n3");
        atomTypes.put(2, "hn");
        atomTypes.put(3, "hn");
        atomTypes.put(4, "hn");
        return atomTypes;
    }

    public Map<Integer, String> makeNOTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "n3");
        atomTypes.put(2, "oh");
        return atomTypes;
    }
    public Map<Integer, String> makeNOTwoTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "ow");
        atomTypes.put(2, "ho");
        atomTypes.put(3, "ho");
        return atomTypes;
    }

    public Map<Integer, String> makeOTwoTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "oh");
        atomTypes.put(2, "oh");
        atomTypes.put(3, "n3");
        return atomTypes;
    }

    public Map<Integer, String> makeSOTwoTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "s2");
        atomTypes.put(2, "o");
        atomTypes.put(3, "o");
        return atomTypes;
    }

    public Map<Integer, String> makeWaterTypes(){
        Map<Integer, String> atomTypes = new HashMap<>();
        atomTypes.put(1, "ow");
        atomTypes.put(2, "ho");
        atomTypes.put(3, "ho");
        return atomTypes;
    }

    public Map<Integer, String> makeSpecies(String confName){
        Map<Integer, String> atomTypes = new HashMap<>();
        switch (confName){
            case "co":
                makeCOTypes();
            case "co2":
                makeCOTwoTypes();
            case "n2":
                makeNTwoTypes();
            case "o2":
                makeOTwoTypes();
            case "nh3":
                makeAmmoniaTypes();
            case "no":
                makeNOTypes();
            case "no2":
                makeNOTwoTypes();
            case "so2":
                makeSOTwoTypes();
            case "water":
                makeWaterTypes();
            case "ch4":
                makeMethaneTypes();
            case "ethane":
                makeEthaneTypes();
            case "ethene":
                makeEtheneTypes();
            case "propane":
                makePropaneTypes();
            case "propene":
                makePropeneTypes();
            case "butane":
                makeButaneTypes();
            case "butadiene":
                makeButadieneTypes();
            case "methanol":
                makeMethanolTypes();
            case "ethanol":
                makeEthanolTypes();
            case "1propanol":
                makeOnePropanolTypes();
            case "2propanol":
                makeTwoPropanolTypes();
            case "isobutanol":
                makeIsoButanolTypes();
            case "benzene":
                makeBenzeneTypes();
            case "toluene":
                makeTolueneTypes();
            case "ethylbenzene":
                makeEthylBenzeneTypes();
            case "oxylene":
                makeOXyleneTypes();
            case "mxylene":
                makeOXyleneTypes();
            case "pxylene":
                makePXyleneTypes();
        }
        return atomTypes;
    }
}
