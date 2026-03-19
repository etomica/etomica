package etomica.GasMOP;

import etomica.space3d.Vector3D;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class PDBWriterNew {


    public static void writePDB(
            String filename,
            Map<Integer, String> atomMap,
            Map<Integer, Vector3D> positionMap,
            ArrayList<ArrayList<Integer>> connectivity
    ) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filename))) {
            // Write HETATM lines
            for (Map.Entry<Integer, String> entry : atomMap.entrySet()) {
                int atomSerial = entry.getKey();
                String element = entry.getValue();
                etomica.space3d.Vector3D pos =  positionMap.get(atomSerial);
                String atomName = element; // Could be adjusted if needed


                //String atomName;
                if (element.length() == 1) {
                    atomName = String.format(" %-2s", element); // ends at column 14
                } else {
                    atomName = String.format("%-3s", element);  // also ends at column 14
                }

                if (element.length() == 1) {
                    atomName = String.format(" %-2s", element); // ends at column 14
                } else {
                    atomName = String.format("%-3s", element);  // ends at column 14
                }


// Right-aligned element symbol at columns 77–78
                String elementField = String.format("%2s", element); // just 2-character string
                if (element != null && !element.isEmpty()) {
                    if (element.length() == 1) {
                        elementField = String.format("%2s", element.toUpperCase());
                    } else {
                        elementField = String.format("%2s", element.substring(0, 1).toUpperCase() + element.substring(1).toLowerCase());
                    }
                } else {
                    elementField = "  "; // Or some other default if the element is missing
                }
               /* String atomLine = String.format(
                        "%-6s%5d %-2s  %3s %1s %4d    %8.3f %8.3f %8.3f%6.2f%6.2f          %2s",
                        "HETATM",               // cols 1–6
                        atomSerial + 1,       // cols 7–11
                        atomName.trim(),        // cols 13–14
                        "UNL",                  // cols 18-20 (shifted right by adding a space before)
                        "",                     // col 21 (blank)
                        atomSerial + 1,       // cols 23-26 (Residue sequence number shifted right by a space)
                        pos.getX(0), pos.getX(1), pos.getX(2),    // cols 31–54
                        1.00,                   // cols 55–60
                        0.00,                   // cols 61–66
                        elementField            // cols 77–78
                );*/
                String atomLine = String.format(
                        "%-6s%5d %-2s   %3s%1s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s",
                        "HETATM",               // cols 1–6
                        atomSerial + 1,       // cols 7–11
                        atomName.trim(),        // cols 13–14
                        "UNL",                  // cols 17-19 (shifted right by the two spaces before)
                        "",                     // col 20
                        atomSerial + 1,       // cols 22-25 (shifted right by the space before)
                        pos.getX(0),           // cols 31-38
                        pos.getX(1),           // cols 39-46
                        pos.getX(2),           // cols 47-54
                        1.00,                   // cols 55–60
                        0.00,                   // cols 61–66
                        elementField            // cols 77–78
                );
                System.out.println(atomLine);
                writer.write(atomLine);
                writer.newLine();
            }

            Map<Integer, Set<Integer>> connectMap = new TreeMap<>();

            for (ArrayList<Integer> connection : connectivity) {
                if (connection.size() < 2) continue;
                int atom1 = connection.get(0);
                for (int i = 1; i < connection.size(); i++) {
                    int atom2 = connection.get(i);
                    // Bidirectional bonding
                    connectMap.computeIfAbsent(atom1, k -> new TreeSet<>()).add(atom2);
                    connectMap.computeIfAbsent(atom2, k -> new TreeSet<>()).add(atom1);
                }
            }

// Write grouped CONECT lines
            for (Map.Entry<Integer, Set<Integer>> entry : connectMap.entrySet()) {
                int atom = entry.getKey();
                Set<Integer> bonded = entry.getValue();
                StringBuilder conectLine = new StringBuilder(String.format("CONECT%5d", atom));
                for (int bondedAtom : bonded) {
                    conectLine.append(String.format("%5d", bondedAtom));
                }
                writer.write(conectLine.toString());
                writer.newLine();
            }

            writer.write("END\n");
        }
    }

    private static String formatAtomName(String name, String element) {
        return (element.length() == 1) ? String.format(" %-3s", name) : String.format("%-4s", name);
    }

    private static String formatElement(String element) {
        return String.format("%2s", element.toUpperCase());
    }
}
