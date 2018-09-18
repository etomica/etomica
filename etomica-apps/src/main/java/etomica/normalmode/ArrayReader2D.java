/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.util.Resources;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;


/** 
 * Reads in lots of 2D arrays from a file.  Returns a 3D array with each n rows of
 * of the file corresponding to one (top-level) element (a 2D array).
 * @author Andrew Schultz
 */
public class ArrayReader2D {
    
    /**
     * Returns 2D arrays from the given file.  The arrays are assumed to be
     * square.  If the given filename is of the form something.bin, it is
     * assume to be in binary format (serialized) and will be deserialized.
     */
    public static double[][][] getFromFile(String fn){
        if (!(fn.matches(".*.bin"))) {
            File binFile = new File(fn+".bin");
            if (binFile.exists()) {
                return getFromBinFile(fn+".bin");
            }
        }
        else if (fn.matches(".*.bin")) {
            return getFromBinFile(fn);
        }
        return getFromFile(fn, 0);
    }

    /**
     * Returns 2D arrays from the given file.  The arrays are rows x columns,
     * where the rows are given and columns is taken as the actual number of
     * columns in the file.
     */
    public static double[][][] getFromFile(String fn, int rows){
        double[][][] S = null;
        try {
            BufferedReader bufReaderS = Resources.openResourceFile(fn, ArrayReader2D.class);
            String line = bufReaderS.readLine();
            ArrayList allS = new ArrayList();
            int dim = 0;
            while (line != null) {
                String[] valueStrings = line.split(" +");
                if (dim == 0) {
                    dim = valueStrings.length;
                    if (rows == 0) {
                        rows = dim;
                    }
                }
                else if (valueStrings.length != dim) {
                    throw new RuntimeException("Not enough data on line "+allS.size()*rows+" of "+fn);
                }
                double[][] thisS = new double[rows][dim];
                for (int i=0; i<rows; i++) {
                    if (i>0) line = bufReaderS.readLine();
                    valueStrings = line.split(" +");
                    for (int j=0; j<dim; j++) {
                        thisS[i][j] = Double.parseDouble(valueStrings[j]);
                    }
                }
                allS.add(thisS);
                line = bufReaderS.readLine();
            }
            S = new double[allS.size()][][];
            for (int i=0; i<S.length; i++) {
                S[i] = ((double[][])allS.get(i));
            }
        }
        catch (IOException e) {
            throw new RuntimeException("Couldn't read data from file "+fn);
        }
        return S;
    }

    /**
     * Deserializes a 3D array from the given file.
     */
    public static double[][][] getFromBinFile(String fn){
        try {
            FileInputStream fis = new FileInputStream(fn);
            ObjectInputStream in = new ObjectInputStream(fis);
            double[][][] bar = (double[][][])in.readObject();
            in.close();
            fis.close();
            return bar;
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }
}
