package etomica.GasMOP;

import java.io.*;
import java.util.*;

public class LammpsToXYZ {

    /*public static void main(String[] args) throws IOException {

        String inputFile = "D:\\Sem-X\\GO\\amber\\test.data";   // your LAMMPS file
        String outputFile = "D:\\Sem-X\\GO\\amber\\test.xyz";

        List<String[]> atoms = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(inputFile))) {

            String line;
            boolean atomsSection = true; // set true if you're directly reading only atom lines

            while ((line = br.readLine()) != null) {

                line = line.trim();

                if (line.isEmpty()) continue;

                // split line
                String[] tokens = line.split("\\s+");

                // Expect at least: id mol type charge x y z
                if (tokens.length < 7) continue;

                String x = tokens[4];
                String y = tokens[5];
                String z = tokens[6];

                atoms.add(new String[]{x, y, z});
            }
        }

        // Write XYZ file
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile))) {

            bw.write(String.valueOf(atoms.size()));
            bw.newLine();
            bw.write("Converted from LAMMPS data file");
            bw.newLine();

            for (String[] atom : atoms) {
                bw.write(String.format("C %s %s %s", atom[0], atom[1], atom[2]));
                bw.newLine();
            }
        }

        System.out.println("XYZ file written: " + outputFile);
    }*/


    /*public static void main(String[] args) throws IOException {

        String inputFile = "D:\\Sem-X\\GO\\amber\\bonds.data";   // your bonds section
        String outputFile = "D:\\Sem-X\\GO\\amber\\bondsOut.txt";

        List<String> validBonds = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(inputFile))) {

            String line;

            while ((line = br.readLine()) != null) {

                line = line.trim();
                if (line.isEmpty()) continue;

                String[] tokens = line.split("\\s+");

                // Expect: id type atom1 atom2
                if (tokens.length < 3) continue;

                int atom1 = Integer.parseInt(tokens[2]);
                int atom2 = Integer.parseInt(tokens[3]);

                // Keep only if both atoms <= 998
                if (atom1 <= 998 && atom2 <= 998) {
                    validBonds.add(line);
                }
                System.out.println(line);
            }
        }

        // Write filtered bonds
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile))) {

            for (String bond : validBonds) {
                bw.write(bond);
                bw.newLine();
            }
        }

        System.out.println("Filtered bonds written: " + outputFile);
        System.out.println("Total valid bonds: " + validBonds.size());
    }*/


    /*public static void main(String[] args) throws IOException {

        String inputFile = "D:\\Sem-X\\GO\\amber\\bonds.data";   // your bonds section
        String outputFile = "D:\\Sem-X\\GO\\amber\\bondsOut.txt";

        List<String[]> bonds = new ArrayList<>();

        // Read bonds
        try (BufferedReader br = new BufferedReader(new FileReader(inputFile))) {
            String line;

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                String[] tokens = line.split("\\s+");

                if (tokens.length < 4) continue;

                bonds.add(tokens);
            }
        }

        // Write fixed bonds
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile))) {

            int newID = 1;

            for (String[] tokens : bonds) {


                int atom1 = Integer.parseInt(tokens[2]);
                int atom2 = Integer.parseInt(tokens[3]);


                // Force:
                // ID → sequential
                // Type → 1
                bw.write(newID + " 1 " + atom1 + " " + atom2);
                bw.newLine();

                newID++;
            }
        }

        System.out.println("Fixed bonds written to: " + outputFile);
    }*/


    public static void main(String[] args) throws IOException {
        String inputFile = "D:\\Sem-X\\GO\\amber\\test.data";   // your bonds section
        String outputFile = "D:\\Sem-X\\GO\\amber\\testOut.txt";

        try (BufferedReader br = new BufferedReader(new FileReader(inputFile));
             BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile))) {

            String line;
            int newID = 1;

            while ((line = br.readLine()) != null) {

                line = line.trim();
                if (line.isEmpty()) continue;

                String[] tokens = line.split("\\s+");

                if (tokens.length < 4) continue;

                double atom1 = Double.parseDouble(tokens[3]);
                double atom2 = Double.parseDouble(tokens[4]);
                double atom3 = Double.parseDouble(tokens[5]);
                double atom4 = Double.parseDouble(tokens[6]);
                // 🔥 Step 1: filter atoms > 998
               // if (atom1 <= 998 && atom2 <= 998 && atom3 <= 998 && atom4 <= 998) {

                    // 🔥 Step 2 + 3:
                    // - renumber ID
                    // - force bond type = 1\
                    bw.write(newID + " 1 " + " 1 "  + atom1+ " " + atom2+ " " + atom3+" " + atom4+" 0 0 0");
                  //  bw.write(newID + " 1 " + atom1 + " " + atom2+ " " + atom3+ " " + atom4);
                    bw.newLine();

                    newID++;
              //  }
            }

            System.out.println("Clean bonds written: " + outputFile);
            System.out.println("Total bonds kept: " + (newID - 1));
        }
    }
}
