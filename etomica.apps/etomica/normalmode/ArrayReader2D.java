package etomica.normalmode;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


/** 
 * Reads in lots of 2D arrays from a file.  Returns a 3D array with each n rows of
 * of the file corresponding to one (top-level) element (a 2D array).
 * @author Andrew Schultz
 */
public class ArrayReader2D {
    
    /**
     * Returns 2D arrays from the given file.  The arrays are assumed to be
     * square.
     */
    public static double[][][] getFromFile(String fn){
        return getFromFile(fn, 0);
    }

    /**
     * Returns 2D arrays from the given file.  The arrays are rows x columns,
     * where the rows are given and columns is taken as the actual number of
     * columns in the file.
     */
    public static double[][][] getFromFile(String fn, int rows){
        FileReader fileReaderS;
        try {
            fileReaderS = new FileReader(fn);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fn+", caught IOException: " + e.getMessage());
        }
        double[][][] S = null;
        try {
            BufferedReader bufReaderS = new BufferedReader(fileReaderS);
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
}
