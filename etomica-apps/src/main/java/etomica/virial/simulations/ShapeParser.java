package etomica.virial.simulations;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import etomica.space.Vector;
import etomica.space.Space;

public class ShapeParser {

    public static ShapeData doParse(String filename, Space space) {
        ShapeData shape = new ShapeData();
        try {
            FileReader fr = new FileReader(filename);
            BufferedReader bufReader = new BufferedReader(fr);
            String line = null;
            while ((line=bufReader.readLine()) != null) {
                line = line.trim();
                String[] spl = line.split(" +");
                if (spl[0].equals("V")) {
                    shape.volume = Double.parseDouble(spl[1]);
                    continue;
                }
                else if (spl[0].equals("vEHS")) {
                    shape.vEHS= Double.parseDouble(spl[1]);
                    continue;
                }
                double[] xyz = new double[3];
                for (int i=0; i<3; i++) {
                    xyz[i] = Double.parseDouble(spl[i]);
                }
                shape.vertices.add(space.makeVector(xyz));
            }
            bufReader.close();
            shape.B2 = 4*shape.volume*shape.vEHS;
            return shape;
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public static class ShapeData {
        public List<Vector> vertices = new ArrayList<Vector>();
        public double volume, vEHS, B2;
    }
}
