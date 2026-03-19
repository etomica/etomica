package etomica.GasMOP;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class radialDist {
    public static void main(String[] args) throws IOException {
        int numBins = 12; // updated to 12x12x12
        double[][][] histogram = new double[numBins][numBins][numBins];

        // Open the file
        BufferedReader br = new BufferedReader(new FileReader("D:\\Sem-VIII\\sorption\\04_27\\methane125M.txt"));

        StringBuilder sb = new StringBuilder();
        String line;
        while ((line = br.readLine()) != null) {
            sb.append(line);
        }
        br.close();

        // Clean up the text
        String data = sb.toString()
                .replace("[", "")   // Remove all [
                .replace("]", "")   // Remove all ]
                .replace(",", " ")  // Replace commas with spaces
                .trim();            // Remove any trailing spaces

        // Now split by spaces
        String[] tokens = data.split("\\s+");

        if (tokens.length != numBins * numBins * numBins) {
            System.out.println("ERROR: Unexpected number of entries!");
            System.out.println("Found: " + tokens.length + ", Expected: " + (numBins*numBins*numBins));
            return;
        }

        // Fill the histogram
        int idx = 0;
        for (int i = 0; i < numBins; i++) {
            for (int j = 0; j < numBins; j++) {
                for (int k = 0; k < numBins; k++) {
                    histogram[i][j][k] = Double.parseDouble(tokens[idx]);
                    idx++;
                }
            }
        }

        System.out.println("Histogram loaded successfully!");
        System.out.println("Sample value: " + histogram[6][6][6]);
    }
}
