package etomica.potential.GAFF;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.potential.UFF.PDBReader;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class PDBReaderOPLS {
    static ISpecies species;
    public static Map<Integer, String> atomMap = new HashMap<>();
    public static HashMap<Integer, String> atomMapModified = new HashMap<>();
    public static ArrayList<ArrayList<Integer>> connectivity = new ArrayList<>();
    public static ArrayList<ArrayList<Integer>> connectivityModified = new ArrayList<>();
    static Map<Integer, String> modifiedAtomIdentifierMapGAFF = new HashMap<>();
    static Map<Integer, String> modifiedAtomIdentifierMapOPLS = new HashMap<>();
    static Map<Integer, Vector> positions = new HashMap<>();
    public static void getSpecies(String confName){
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        SpeciesBuilder speciesBuilderNewMod =  new SpeciesBuilder(Space3D.getInstance());
        PDBReader pdbReader = new PDBReader();
        pdbReader.readPDBFile(confName);
        positions = pdbReader.getPositions();
        atomMap = pdbReader.getAtomMap();
        atomMapModified = pdbReader.getatomMapModified(atomMap);
        connectivity = pdbReader.getConnectivity();
        connectivityModified = pdbReader.getconnectivityModified(connectivity);
        System.out.println("\n New One ");
        System.out.println(positions);
        System.out.println(atomMap);
        System.out.println(atomMapModified);
        System.out.println(connectivity);
        System.out.println(connectivityModified);
        getModifiedAtomIdentifierMapGAFF(connectivityModified, atomMapModified);
    }

    public static Map<Integer,String> getModifiedAtomIdentifierMapGAFF (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified) {
        Map<Integer, AtomType> atomIdentifierMapGAFF = atomIdentifierGAFF(connectivityModified, atomMapModified);
        Map<Integer, AtomType> atomIdentifierMapOPLS = atomIdentifierOPLS(connectivityModified, atomMapModified);
        int i =0;
        for (Map.Entry<Integer, AtomType> entry : atomIdentifierMapGAFF.entrySet()) {
            String value = entry.getValue().toString();
           // String valueOPLS = atomIdentifierMapOPLS.get(i).toString();

            value = value.replace("AtomType[", "").replace("]", "");
           // valueOPLS = valueOPLS.replace("AtomType[", "").replace("]", "");
            modifiedAtomIdentifierMapGAFF.put(entry.getKey(), value);
            //modifiedAtomIdentifierMapOPLS.put(i, valueOPLS);
        }
        System.out.println(modifiedAtomIdentifierMapGAFF + " modified Identifier");
        System.out.println(modifiedAtomIdentifierMapOPLS + " modified OPLS");
        return modifiedAtomIdentifierMapGAFF;
    }

    public static Map<Integer,String> getModifiedAtomIdentifierMapOPLS(){
        return modifiedAtomIdentifierMapOPLS;
    }

    protected static  Map<Integer, AtomType> atomIdentifierGAFF(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        Map<Integer, AtomType> atomIdentifierMapOPLS = new HashMap<>();
        int numHydr=0, numCarb=0, numOxy=0;
        System.out.println("GAFF");
        for (int i=0; i<connectivityModified.size(); i++){
            Integer retriveArrayFirstElementNumber = connectivityModified.get(i).get(0);
            String retriveArrayFirstElementName = atomMapModified.get(retriveArrayFirstElementNumber);
            System.out.println(retriveArrayFirstElementName+ " " + retriveArrayFirstElementNumber);
            int arraySize = connectivityModified.get(i).size();
            if(retriveArrayFirstElementName.equals("C")){
                int arrayListSize = connectivityModified.get(i).size();
                for(int j=1; j<=arrayListSize; j++){
                    String connectedElement = atomMapModified.get(connectivityModified.get(i).get(j));
                    if(connectedElement.equals("C")){
                        numCarb++;
                    } else if (connectedElement.equals("H")) {
                        numHydr++;
                    }else {
                        numOxy++;
                    }
                }
                if(arraySize == 5){
                    //single bond carbon
                    if(numCarb ==0){
                        AtomType C4ane = new AtomType(Carbon.INSTANCE, "C4ane");
                        atomIdentifierMapOPLS.put(i, C4ane);
                    } else if (numCarb ==1) {
                        AtomType C3ane = new AtomType(Carbon.INSTANCE, "C3ane");
                        atomIdentifierMapOPLS.put(i, C3ane);
                    } else if (numCarb ==2) {
                        AtomType C2ane = new AtomType(Carbon.INSTANCE, "C2ane");
                        atomIdentifierMapOPLS.put(i, C2ane);
                    } else if (numCarb ==3) {
                        AtomType C1ane = new AtomType(Carbon.INSTANCE, "C1ane");
                        atomIdentifierMapOPLS.put(i, C1ane);
                    }
                } else if (arraySize == 4) {
                    //double bond Carbon
                    String retriveArrayOtherElementName;
                    for(int j =0; j<arraySize; j++){
                        retriveArrayOtherElementName = atomMapModified.get(connectivityModified.get(i).get(j));
                        if(retriveArrayFirstElementName.equals("O")){
                            AtomType c = AtomType.element(Carbon.INSTANCE, "c");
                            atomIdentifierMapOPLS.put(i, c);
                        } else if (retriveArrayOtherElementName.equals("C")) {
                            AtomType c2 = AtomType.element(Carbon.INSTANCE, "c2");
                            atomIdentifierMapOPLS.put(i, c2);
                        }
                    }
                } else if (arraySize == 3) {
                    //triple bonds
                    AtomType c1 = AtomType.element(Carbon.INSTANCE, "c1");
                    atomIdentifierMapOPLS.put(i, c1);
                } else {
                    //for CO or CS
                    AtomType c = AtomType.element(Carbon.INSTANCE, "c");
                    atomIdentifierMapOPLS.put(i, c);
                }

            } else if (retriveArrayFirstElementName.equals("H")) {
                
            } else if (retriveArrayFirstElementName.equals("O")) {
                
            } else if (retriveArrayFirstElementName.equals("N")) {
                
            }
        }
        return atomIdentifierMapOPLS;
    }
    protected static  Map<Integer, AtomType> atomIdentifierOPLS(ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer, String> atomMapModified){
        Map<Integer, AtomType> atomIdentifierMapOPLS = new HashMap<>();
       /* System.out.println("OPLS");
        for (int i=0; i<connectivityModified.size(); i++){
            System.out.println(connectivityModified.get(i));
        }*/
        return atomIdentifierMapOPLS;
    }
    public static double [] nonBondedPot(String atomtype){
        HashMap<String, double[]> atomicConstant = new HashMap<>();
        atomicConstant.put("C4ane", new double[]{-0.240,3.5,0.066});
        atomicConstant.put("C3ane", new double[]{ -0.180,3.5,0.066 });
        atomicConstant.put("C2ane", new double[]{-0.0120 ,3.5 ,0.066 });
        atomicConstant.put("C1ane", new double[]{ -0.060,3.5 ,0.066 });
        atomicConstant.put("Hane", new double[]{ 0.06,2.5 ,0.030 });
        atomicConstant.put("Carom", new double[]{ -0.115,3.55 ,0.070 });
        atomicConstant.put("Harom", new double[]{ 0.115,2.42 ,0.03 });
        atomicConstant.put("Ctol", new double[]{ -0.065,3.5,0.066 });
        atomicConstant.put("Cethylbenz", new double[]{-0.005,3.5,0.066 });
        atomicConstant.put("C0ene", new double[]{ 0,3.55,0.076 });
        atomicConstant.put("C1ene", new double[]{ -0.115,3.55 ,0.076 });
        atomicConstant.put("C2ene", new double[]{ -0.230,3.55 ,0.076 });
        atomicConstant.put("Hene", new double[]{ 0.115,2.42 ,0.030 });
        atomicConstant.put("Oalc", new double[]{ -0.683,3.12 ,0.17 });
        atomicConstant.put("Halc", new double[]{ 0.418, 0 ,0 });
        atomicConstant.put("HalcC", new double[]{ 0.040,2.5 ,0.030 });
        atomicConstant.put("C32alc", new double[]{ 0.145,3.5 ,0.066 });
        atomicConstant.put("C1alc", new double[]{ 0.205, 3.5,0.066 });
        atomicConstant.put("C0alc", new double[]{ 0.265, 3.5,0.066 });
        double [] sample = atomicConstant.get(atomtype);
        return sample;
    }


    public static void main(String[] args) {
        String confName = "F://Avagadro//molecule//ethane";
        getSpecies(confName);
    }
}
