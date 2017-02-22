/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;


/**
 * Reads a file containing n mxm matrices as (serialized) binary data and
 * writes them to another file as text (ASCII).
 * 
 * @author Andrew Schultz
 */
public class ArrayB2A {

    public static void main(String[] args) {

        
        if (args.length < 1) {
            throw new RuntimeException("usage: B2A filename");
        }
        double[][][] x;
        try {
            FileInputStream fis = new FileInputStream(args[0]);
            ObjectInputStream in = new ObjectInputStream(fis);
            x = (double[][][])in.readObject();
            in.close();
            fis.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        }

        try {
            FileWriter fw = new FileWriter(args[0]+".txt");
            for (int i = 0; i<x.length; i++) {
                for (int j = 0; j < x[i].length; j++) {
                    fw.write(Double.toString(x[i][j][0]));
                    for (int k = 1; k < x[i][j].length; k++) {
                        fw.write(" " + x[i][j][k]);
                    }
                    fw.write("\n");
                }
            }
            fw.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }

    }
}
