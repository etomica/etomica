/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.space.Space;
import etomica.space.Vector;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

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
                shape.vertices.add(Vector.of(xyz));
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
