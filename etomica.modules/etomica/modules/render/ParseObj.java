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
        ArrayList<Vector3D> vertices = new ArrayList();
        try {
            fileReader = new FileReader(fileName);
        } catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            for (int iLeaf=0; iLeaf<30; iLeaf++) {
                String currentLine = bufReader.readLine();
                System.out.println(currentLine);
//                vertices.add(vertex);
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        positions = (Vector3D[])vertices.toArray();

    }
    
    public int nAtoms;
    public Vector3D[] positions;
    
    public static void main(String[] args) {
        String file = "/Users/kofke/Documents/workspace/car self-assembly/mustang.txt";
        ParseObj parser = new ParseObj(file);

    }

}
