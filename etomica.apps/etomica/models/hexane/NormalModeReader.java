package etomica.models.hexane;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


/** 
 * Reads in a file.
 * @author cribbin
 *
 */
public class NormalModeReader {
    
    public void readFromFile(String fn){
        FileReader fileReaderS;
        try {
            fileReaderS = new FileReader(fn+".S");
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fn+", caught IOException: " + e.getMessage());
        }
        S = null;
        try {
            BufferedReader bufReaderS = new BufferedReader(fileReaderS);
            String line = bufReaderS.readLine();
            ArrayList allS = new ArrayList();
            int dim = 0;
            while (line != null) {
                String[] valueStrings = line.split(" +");
                if (dim == 0) {
                    dim = valueStrings.length;
                }
                else if (valueStrings.length != dim) {
                    throw new RuntimeException("Not enough data on line "+allS.size()*dim+" of "+fn);
                }
                double[][] thisS = new double[dim][dim];
                for (int i=0; i<dim; i++) {
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
                for (int j=0; j<dim; j++) {
                    for (int k=0; k<dim; k++) {
                        S[i] = ((double[][])allS.get(i));
                    }
                }
            }
        }
        catch (IOException e) {
            throw new RuntimeException("Couldn't read data from file "+fn);
        }
    }
    
    public double[][][] S;
    
    public static void main(String[] args) {
        NormalModeReader reader = new NormalModeReader();
        reader.readFromFile(args[0]);
        double[][][] S = reader.S;
        for (int i=0; i<S.length; i++) {
            
            for (int j=0; j<S[i].length; j++) {
                for (int k=0; k<S[i][j].length; k++) {
                    System.out.print(S[i][j][k]+" ");
                }
                System.out.println("");
            }
            System.out.println("");
        }
    }
}
