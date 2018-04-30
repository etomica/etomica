/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import etomica.space.Vector;
import etomica.space.Space;


/** 
 * Reads in lots of 1D arrays from a file.  Returns a 2D array with each row of
 * of the file corresponding to one (top-level) element.  When called via 
 * getVectorsFromFile, it returns a 1D array of Vectors with each row in the 
 * file corresponding to a Vector.
 * @author Andrew Schultz
 */
public class ArrayReader1D {
    
    public static Vector[] getVectorsFromFile(String fn) {
        return (Vector[])getFromFile(fn, true);
    }
    
    public static double[][] getFromFile(String fn){
        return (double[][])getFromFile(fn, false);
    }
    
    protected static Object[] getFromFile(String fn, boolean useVectors){
        FileReader fileReaderQ;
        try {
            fileReaderQ = new FileReader(fn);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fn+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReaderQ = new BufferedReader(fileReaderQ);
            String line = bufReaderQ.readLine();
            ArrayList allQ = new ArrayList();
            int dim = 0;
            while (line != null) {
                String[] valueStrings = line.split(" +");
                if (dim == 0) {
                    dim = valueStrings.length;
                }
                else if (valueStrings.length != dim) {
                    throw new RuntimeException("Not enough data on line "+allQ.size()+" of "+fn);
                }
                double[] thisQ = new double[dim];
                for (int i=0; i<dim; i++) {
                    thisQ[i] = Double.parseDouble(valueStrings[i]);
                }
                if (useVectors) {
                	Space space = Space.getInstance(thisQ.length);
                    allQ.add(Vector.of(thisQ));
                }
                else {
                    allQ.add(thisQ);
                }
                line = bufReaderQ.readLine();
            }
            
            Object[] q = null;
            if (useVectors) {
                q = new Vector[allQ.size()];
            }
            else {
                q = new double[allQ.size()][];
            }
            for (int i=0; i<q.length; i++) {
                for (int j=0; j<dim; j++) {
                    for (int k=0; k<dim; k++) {
                        q[i] = allQ.get(i);
                    }
                }
            }
            return q;
        }
        catch (IOException e) {
            throw new RuntimeException("Couldn't read data from file "+fn);
        }
    }
}
