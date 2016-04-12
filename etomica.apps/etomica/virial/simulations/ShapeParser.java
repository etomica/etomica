package etomica.virial.simulations;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import etomica.api.IVector;
import etomica.space.ISpace;

public class ShapeParser {

    public static List<IVector> getVertices(String filename, ISpace space) {
        List<IVector> list = new ArrayList<IVector>();
        try {
            FileReader fr = new FileReader(filename);
            BufferedReader bufReader = new BufferedReader(fr);
            String line = null;
            while ((line=bufReader.readLine()) != null) {
                line = line.trim();
                String[] spl = line.split(" +");
                double[] xyz = new double[3];
                for (int i=0; i<3; i++) {
                    xyz[i] = Double.parseDouble(spl[i]);
                }
                list.add(space.makeVector(xyz));
            }
            bufReader.close();
            return list;
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

}
