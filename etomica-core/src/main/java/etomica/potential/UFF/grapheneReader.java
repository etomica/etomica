package etomica.potential.UFF;

import com.fasterxml.jackson.databind.node.ObjectNode;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.chem.elements.*;
import etomica.molecule.IMolecule;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class grapheneReader {
    static ISpecies species;
    double massSum;
    public Map<Integer, String> atomMapGraphene = new HashMap<>();
    public Map<Integer,Vector> positions = new HashMap<>();
    public ArrayList<ArrayList<Integer>> connectivityGrapehene = new ArrayList<>();
    public Map<Integer, String> modifiedAtomIdentifierMapGraphene = new HashMap<>();
    public Map<Integer, Vector> positionsGraphene = new HashMap<>();
    public Map<String, AtomType> elementReceiverMap = new HashMap<>();
    public Map<String, AtomType> typeMap = new HashMap<>();
    public List<int[]> dupletsSorted, tripletsSorted, quadrupletsSorted = new ArrayList<>();
    Map<String, AtomType> typeMapNew = new HashMap<>();

    public ISpecies getSpecies(String confName,Vector vectorCOM, boolean isInfinite){
        AtomType typeNew;
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        SpeciesBuilder speciesBuilderNewMod =  new SpeciesBuilder(Space3D.getInstance());
        PDBReader pdbReader = new PDBReader();
        pdbReader.readPDBFile(confName);
        positionsGraphene = pdbReader.getPositions();
        atomMapGraphene = pdbReader.getAtomMap();
        connectivityGrapehene = pdbReader.getConnectivityForGraphene();
      //  System.out.println(positionsGraphene);
       // System.out.println(atomMapGraphene);
       // System.out.println("\n New One ");
        modifiedAtomIdentifierMapGraphene= getModifiedAtomIdentifierMap(atomMapGraphene);
        setBondsGrapehene(connectivityGrapehene, modifiedAtomIdentifierMapGraphene, dupletsSorted, tripletsSorted, quadrupletsSorted);
        System.exit(1);
        Space space = Space3D.getInstance();
        Vector center = space.makeVector();
        Vector dr = Vector.d(center.getD());
      //  System.out.println(modifiedAtomIdentifierMapGraphene);
        int n =0;
        if(modifiedAtomIdentifierMapGraphene.get(0) == null){
            n =1;
        }
        for(int i = n; i < modifiedAtomIdentifierMapGraphene.size(); i++) {
            String symbol = String.valueOf(modifiedAtomIdentifierMapGraphene.get(i));
            //System.out.println(symbol);
            // System.out.println(atomIdentifierMap.get(i));
            AtomType newName = returnElement(symbol, isInfinite);
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
    public Map<Integer,String> getModifiedAtomIdentifierMap ( Map<Integer, String> atomMapModified) {
        HashSet<String> uniquesAtoms = new HashSet<>();
        for (String str : atomMapModified.values()) {
            uniquesAtoms.add(str);
        }
      //  System.out.println(uniquesAtoms);
        int n=0;
        if(atomMapModified.get(0) ==null){
            n=1;
        }
        for (String str : atomMapModified.values()) {
            n++;
            if(str == null){
                throw new RuntimeException("Does not exist " + n);
            }
        }
        Map<Integer, AtomType> atomIdentifierMapGAFF = atomIdentifierGraphene(atomMapModified);
        int i =0;
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

    protected Map<Integer, AtomType> atomIdentifierGraphene(Map<Integer, String> atomMapModified){
        Map<Integer, AtomType> atomIdentifierMapGraphene = new HashMap<>();
        int n =0;
        if(atomMapModified.get(0) == null){
            n=1;
        }
        int numHydr=0, numCarb=0, numOxy=0;
        for ( int i = n; i<atomMapModified.size(); i++){
            String atomName = atomMapModified.get(i);
            if(atomName.equals("CX")){
                AtomType CX = new AtomType(Carbon.INSTANCE, "CX");
                atomIdentifierMapGraphene.put(i, CX);
            } else if (atomName.equals("CY")) {
                AtomType CY = new AtomType(Carbon.INSTANCE,"CY");
                atomIdentifierMapGraphene.put(i, CY);
            } else if (atomName.equals("HK")) {
                AtomType HK = new AtomType(Hydrogen.INSTANCE,"HK");
                atomIdentifierMapGraphene.put(i, HK);
            } else if (atomName.equals("C4")) {
                AtomType C4 = new AtomType(Carbon.INSTANCE,"C4");
                atomIdentifierMapGraphene.put(i, C4);
            }else if (atomName.equals("OJ")) {
                AtomType OJ = new AtomType(Oxygen.INSTANCE,"OJ");
                atomIdentifierMapGraphene.put(i, OJ);
            }else if (atomName.equals("OK")) {
                AtomType OK = new AtomType(Oxygen.INSTANCE,"OK");
                atomIdentifierMapGraphene.put(i, OK);
            }else if (atomName.equals("OL")) {
                AtomType OL = new AtomType(Oxygen.INSTANCE,"OL");
                atomIdentifierMapGraphene.put(i, OL);
            }else if (atomName.equals("OE")) {
                AtomType OE = new AtomType(Oxygen.INSTANCE,"OE");
                atomIdentifierMapGraphene.put(i, OE);
            }else if (atomName.equals("H1")) {
                AtomType H1 = new AtomType(Hydrogen.INSTANCE,"H1");
                atomIdentifierMapGraphene.put(i, H1);
            }else if (atomName.equals("H13")) {
                AtomType H13 = new AtomType(Hydrogen.INSTANCE,"H13");
                atomIdentifierMapGraphene.put(i, H13);
            }else if (atomName.equals("H15")) {
                AtomType H15 = new AtomType(Hydrogen.INSTANCE,"H15");
                atomIdentifierMapGraphene.put(i, H15);
            }else if (atomName.equals("H5")) {
                AtomType H5 = new AtomType(Hydrogen.INSTANCE,"H5");
                atomIdentifierMapGraphene.put(i, H5);
            }else if (atomName.equals("H6")) {
                AtomType H6 = new AtomType(Hydrogen.INSTANCE,"H6");
                atomIdentifierMapGraphene.put(i, H6);
            }else if (atomName.equals("H7")) {
                AtomType H7 = new AtomType(Hydrogen.INSTANCE,"H7");
                atomIdentifierMapGraphene.put(i, H7);
            }else if (atomName.equals("H8")) {
                AtomType H8 = new AtomType(Hydrogen.INSTANCE,"H8");
                atomIdentifierMapGraphene.put(i, H8);
            }else if (atomName.equals("H9")) {
                AtomType H9 = new AtomType(Hydrogen.INSTANCE,"H9");
                atomIdentifierMapGraphene.put(i, H9);
            }else if (atomName.equals("H10")) {
                AtomType H10 = new AtomType(Hydrogen.INSTANCE,"H10");
                atomIdentifierMapGraphene.put(i, H10);
            }else if (atomName.equals("H11")) {
                AtomType H11 = new AtomType(Hydrogen.INSTANCE,"H11");
                atomIdentifierMapGraphene.put(i, H11);
            }else if (atomName.equals("H12")) {
                AtomType H12 = new AtomType(Hydrogen.INSTANCE,"H12");
                atomIdentifierMapGraphene.put(i, H12);
            }else if (atomName.equals("H14")) {
                AtomType H14 = new AtomType(Hydrogen.INSTANCE,"H14");
                atomIdentifierMapGraphene.put(i, H14);
            }else if (atomName.equals("H2")) {
                AtomType H2 = new AtomType(Hydrogen.INSTANCE,"H2");
                atomIdentifierMapGraphene.put(i, H2);
            }else if (atomName.equals("H3")) {
                AtomType H3 = new AtomType(Hydrogen.INSTANCE,"H3");
                atomIdentifierMapGraphene.put(i, H3);
            }else if (atomName.equals("H4")) {
                AtomType H4 = new AtomType(Hydrogen.INSTANCE,"H4");
                atomIdentifierMapGraphene.put(i, H4);
            }

        }
        return atomIdentifierMapGraphene;
    }

    public void setBondsGrapehene (ArrayList<ArrayList<Integer>> connectivityModified, Map<Integer,String> atomIdentifierModified, List<int[]> dupletsSorted, List<int[]> tripletesSorted, List<int[]> quadrupletsSorted ){
        System.out.println(connectivityModified);
        System.out.println(atomIdentifierModified);
    }


    /*protected static ArrayList<ArrayList<Integer>> getConnectivity ( Map<Integer, Vector> positionsGraphene){
        Map<Integer, ArrayList<Integer>> linkerMap = new HashMap<>();
        List<Integer> connect = new ArrayList<>();
        ArrayList<Integer> temp = new ArrayList<>();
        double targetDifferenceX = 1.228;
        double targetDifferenceY = 0.709;
        for (int i =3; i<positionsGraphene.size(); i++){
            double[] v1 = new double[3];
            positionsGraphene.get(i).assignTo(v1);
            //System.out.println(positionsGraphene.get(i).toArray() [0] + " " +positionsGraphene.get(i).toArray() [1]+ " grapo");
            double value = 0;
            for (int j=i+1; j< positionsGraphene.size(); j++){
                double[] v2 = new double[3];
                positionsGraphene.get(j).assignTo(v2);
                double diffX = Math.abs(v2[0] - v1[0]) - targetDifferenceX;
                double diffY = Math.abs(v2[1] - v1[1]) - targetDifferenceY;
                System.out.println(Arrays.toString(v1) + " " + Arrays.toString(v2) + " " + diffX + " " + diffY);
                if(diffX < 1E-6 && diffY < 1E-6){
                    connect.add(i);
                    connect.add(j);
                    System.out.println(connect);
                    connect.clear();
                }
            }
        }
        return connectivityGrapehene;
    }*/
    public AtomType returnElement(String elementName, boolean isInfinite){
        if(isInfinite){
            switch (elementName){
                case "CX":
                    AtomType CX = new AtomType(new ElementSimple("CX", Double.POSITIVE_INFINITY), "CX");
                    elementReceiverMap.put("CX", CX);
                case "CY":
                    AtomType CY = new AtomType(new ElementSimple("CY", Double.POSITIVE_INFINITY), "CY");
                    elementReceiverMap.put("CY", CY);
                case "C4":
                    AtomType C4 = new AtomType(new ElementSimple("C4", Double.POSITIVE_INFINITY), "C4");
                    elementReceiverMap.put("C4", C4);
                case "HK":
                    AtomType HK = new AtomType(new ElementSimple("HK", Double.POSITIVE_INFINITY), "HK");
                    elementReceiverMap.put("HK", HK);
            }

        } else {
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
        }
        return elementReceiverMap.get(elementName);
    }
    public double[] atomicPotGraphene (String atomType){
        HashMap<String, double[]> atomicConstant = new HashMap<>();
        atomicConstant.put("CA", new double[]{0,0, 3.851,  0.105});
        atomicConstant.put("CY", new double[]{0,0,3.851 , 0.105});
        atomicConstant.put("CX", new double[]{3.42, 0.170});//1.9080 , 0.0860
        atomicConstant.put("C4", new double[]{3.55, 0.07});
        atomicConstant.put("C6", new double[]{3.851, 0.105});
        atomicConstant.put("C5", new double[]{3.851, 0.105});
        atomicConstant.put("C2", new double[]{3.851, 0.105});
        atomicConstant.put("CZ", new double[]{3.851, 0.105});
        atomicConstant.put("C1", new double[]{3.851, 0.105});
        atomicConstant.put("C3", new double[]{3.851, 0.105});
        atomicConstant.put("C", new double[]{3.5 , 0.105 });
        atomicConstant.put("OE", new double[]{3.5, 0.2104});
        atomicConstant.put("O2", new double[]{3.5, 0.2104});
        atomicConstant.put("OJ", new double[]{ 3.5, 0.2104});
        atomicConstant.put("OK", new double[]{ 3.5, 0.2104 });
        atomicConstant.put("OL", new double[]{ 3.5, 0.2104 });
        atomicConstant.put("OS", new double[]{3.5,  0.1700 });
        atomicConstant.put("OH", new double[]{ 3.5,  0.2104});
        atomicConstant.put("O", new double[]{3.5 , 0.2100  });
        atomicConstant.put("O1", new double[]{3.5 , 0.2100  });
        atomicConstant.put("HK", new double[]{2.886, 0.044 });
        atomicConstant.put("H1", new double[]{2.886, 0.044 });
        atomicConstant.put("H2", new double[]{2.886, 0.044 });
        atomicConstant.put("H3", new double[]{2.886, 0.044 });
        atomicConstant.put("H4", new double[]{2.886, 0.044 });
        atomicConstant.put("H5", new double[]{2.886, 0.044 });
        atomicConstant.put("H6", new double[]{2.886, 0.044 });
        atomicConstant.put("H7", new double[]{2.886, 0.044 });
        atomicConstant.put("H8", new double[]{2.886, 0.044 });
        atomicConstant.put("H9", new double[]{2.886, 0.044 });
        atomicConstant.put("H10", new double[]{2.886, 0.044 });
        atomicConstant.put("H11", new double[]{2.886, 0.044 });
        atomicConstant.put("H12", new double[]{2.886, 0.044 });
        atomicConstant.put("H13", new double[]{2.886, 0.044 });
        atomicConstant.put("H14", new double[]{2.886, 0.044 });
        atomicConstant.put("H15", new double[]{2.886, 0.044 });
        atomicConstant.put("N2", new double[]{ 3.66, 0.069});
        atomicConstant.put("N1", new double[]{ 3.66, 0.069});
        atomicConstant.put("N3", new double[]{ 3.66, 0.069});
        double [] sample = atomicConstant.get(atomType);
        return sample;
    }

    public static void main(String[] args) {
        grapheneReader grapheneReader = new grapheneReader();
        String confName = "F://Avagadro//GO//pristine//pris55";
        double[] arrayCOM= new double[]{10,10,10};
        Vector vectorCOM = Vector.of(arrayCOM);
        ISpecies graphene = grapheneReader.getSpecies(confName, vectorCOM, false);
    }
}
