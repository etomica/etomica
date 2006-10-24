package etomica.eigenstuff;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


/** 
 * Reads in a file.
 * @author cribbin
 *
 */
public class ReaderIn {
    
    public static double[][][] getFromFile(String fn, int dim){
        FileReader fileReader;
        try {
            fileReader = new FileReader(fn);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fn+", caught IOException: " + e.getMessage());
        }
        double[][][] valuesFromFile = null;
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            String line = bufReader.readLine();
            ArrayList allValues = new ArrayList();
            while (line != null) {
                String[] valueStrings = line.split(" +");
                if (valueStrings.length != dim*dim) {
                    throw new RuntimeException("Not enough data on line "+allValues.size()+" of "+fn);
                }
                double[][] values = new double[dim][dim];
                for (int i=0; i<dim; i++) {
                    for (int j=0; j<dim; j++) {
                        values[i][j] = Double.parseDouble(valueStrings[i*dim+j]);
                    }
                }
                allValues.add(values);
                line = bufReader.readLine();
            }
            valuesFromFile = new double[allValues.size()][][];
            for (int i=0; i<valuesFromFile.length; i++) {
                for (int j=0; j<dim; j++) {
                    for (int k=0; k<dim; k++) {
                        valuesFromFile[i] = ((double[][])allValues.get(i));
                    }
                }
            }
        }
        catch (IOException e) {
            throw new RuntimeException("Couldn't read data from file "+fn);
        }
        return valuesFromFile;
    }
    
    public static void main(String[] args) {
        double[][][] values = getFromFile(args[0], 3);
        for (int i=0; i<values.length; i++) {
            for (int j=0; j<values[i].length; j++) {
                for (int k=0; k<values[i][j].length; k++) {
                    System.out.print(values[i][j][k]+" ");
                }
            }
            System.out.println("");
        }
    }
}
