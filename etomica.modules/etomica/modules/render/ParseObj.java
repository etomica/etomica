package etomica.modules.render;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import etomica.api.IAtom;
import etomica.space3d.Vector3D;

public class ParseObj {
    
    public ParseObj(String fileName) {
        
        FileReader fileReader;
        
        int countV = 0;
        try {
            fileReader = new FileReader(fileName);
        } catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            String currentLine;
//            for (int iLeaf=0; iLeaf<30; iLeaf++) {
            int[] indexRef = new int[100];
            double tol = 1e-14;
            while((currentLine = bufReader.readLine()) != null) {
                String[] coordStr = currentLine.split(" ");
                if(coordStr[0].equals("v") /*&& countV < 1000*/) {
                    Vector3D vertex = new Vector3D(Double.valueOf(coordStr[1]).doubleValue(),Double.valueOf(coordStr[2]).doubleValue(),Double.valueOf(coordStr[3]).doubleValue());
                    if(countV == indexRef.length) indexRef = expand(indexRef);
                    
                    //look for vertex position in list of ones read previously
                    boolean foundVertex = false;
                    for(int i=0; i<vertices.size(); i++) {
                        double r2 = vertex.Mv1Squared(vertices.get(i));
                        if(r2 < tol) {//found vertex in list
                            foundVertex = true;
                            indexRef[countV] = i;
                            break;
                        }
                    }
                    if(!foundVertex) {
                        indexRef[countV] = vertices.size();
                        vertices.add(vertex);
                    }
                    countV++;
                }
                if(coordStr[0].equals("vn")) {
                    //vertex normal; do not use
                }
                if(coordStr[0].equals("f")) {
                    try {
                    int i0 = indexRef[Integer.valueOf(coordStr[1].split("//")[0]).intValue() - 1];
                    int i1 = indexRef[Integer.valueOf(coordStr[2].split("//")[0]).intValue() - 1];
                    int i2 = indexRef[Integer.valueOf(coordStr[3].split("//")[0]).intValue() - 1];
                    Vector3D v0 = vertices.get(i0);                    
                    Vector3D v1 = vertices.get(i1);                    
                    Vector3D v2 = vertices.get(i2);                    
                    bondList.add(new BondInfo(i0, i1, v0.Mv1Squared(v1)));
                    bondList.add(new BondInfo(i0, i2, v0.Mv1Squared(v2)));
                    bondList.add(new BondInfo(i1, i2, v1.Mv1Squared(v2)));
                    } catch(ArrayIndexOutOfBoundsException ex) {};
                }
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        nAtoms = vertices.size();
        System.out.println(nAtoms+" atoms from "+countV+" vertices");

    }
    
    private int[] expand(int[] index) {
        int[] newIndex = new int[index.length+100];
        System.arraycopy(index, 0, newIndex, 0, index.length);
        return newIndex;
    }
    
    public int nAtoms;
    
    public ArrayList<BondInfo> bondList = new ArrayList<BondInfo>();
    ArrayList<Vector3D> vertices = new ArrayList<Vector3D>();
    
    public class BondInfo {
        public BondInfo(int i0, int i1, double bondLengthSquared) {
            this.i0 = i0;
            this.i1 = i1;
            this.bondLengthSquared = bondLengthSquared;
        }
        int i0, i1;
        double bondLengthSquared;
    }
    public static void main(String[] args) {
        String file = "/Users/kofke/Documents/workspace/car self-assembly/mustang.txt";
        ParseObj parser = new ParseObj(file);

    }

}
