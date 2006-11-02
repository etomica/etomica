package etomica.models.hexane;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


/** 
 * Reads in the eignavlues from a file and returns omega/(kT)^(1/2).
 * @author Andrew Schultz
 */
public class OmegaReader {
    
    public static double[][] getFromFile(String fn){
        String filename = fn+".val";
        FileReader fileReaderVal;
        try {
            fileReaderVal = new FileReader(filename);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+filename+", caught IOException: " + e.getMessage());
        }
        double[][] omega = null;
        try {
            BufferedReader bufReaderS = new BufferedReader(fileReaderVal);
            String line = bufReaderS.readLine();
            ArrayList allV = new ArrayList();
            int dim = 0;
            while (line != null) {
                String[] valueStrings = line.split(" +");
                if (dim == 0) {
                    dim = valueStrings.length;
                }
                else if (valueStrings.length != dim) {
                    throw new RuntimeException("Not enough data on line "+allV.size()+" of "+filename);
                }
                double[] thisV = new double[dim];
                valueStrings = line.split(" +");
                for (int j=0; j<dim; j++) {
                    thisV[j] = Double.parseDouble(valueStrings[j]);
                }
                allV.add(thisV);
                line = bufReaderS.readLine();
            }
            omega = new double[allV.size()][dim];
            for (int i=0; i<omega.length; i++) {
                for (int j=0; j<dim; j++) {
                   omega[i][j] = 1/((double[])allV.get(i))[j];
                }
            }
        }
        catch (IOException e) {
            throw new RuntimeException("Couldn't read data from file "+filename);
        }
        return omega;
    }
    
    public static void main(String[] args) {
        String filename = "normal_modes";
        if (args.length > 0) {
            filename = args[0];
        }
        double[][] omega = OmegaReader.getFromFile(filename);
        double f = 0;
        for (int i=0; i<omega.length; i++) {
            for (int j=0; j<omega[i].length; j++) {
                f += Math.log(omega[i][j]);
            }
        }
        //XXX appropriate constants need to be used here
        System.out.println("f="+f);
    }
}
