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
        FileReader fileReader;
        try {
            fileReader = new FileReader(fn);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fn+", caught IOException: " + e.getMessage());
        }
        S = null;
        q = null;
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            String line = bufReader.readLine();
            ArrayList allQ = new ArrayList();
            ArrayList allS = new ArrayList();
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
                for (int j=0; j<dim; j++) {
                    thisQ[j] = Double.parseDouble(valueStrings[j]);
                }
                double[][] thisS = new double[dim][dim];
                for (int i=0; i<dim; i++) {
                    line = bufReader.readLine();
                    valueStrings = line.split(" +");
                    for (int j=0; j<dim; j++) {
                        thisS[i][j] = Double.parseDouble(valueStrings[j]);
                    }
                }
                allQ.add(thisQ);
                allS.add(thisS);
                line = bufReader.readLine();
            }
            q = new double[allQ.size()][];
            S = new double[allQ.size()][][];
            for (int i=0; i<q.length; i++) {
                q[i] = (double[])allQ.get(i);
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
    
    public double[][] q;
    public double[][][] S;
    
    public static void main(String[] args) {
        NormalModeReader reader = new NormalModeReader();
        reader.readFromFile(args[0]);
        double[][] q = reader.q;
        double[][][] S = reader.S;
        for (int i=0; i<q.length; i++) {
            for (int k=0; k<q[i].length; k++) {
                System.out.print(q[i][k]+" ");
            }
            System.out.println("");
            
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
