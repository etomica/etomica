package etomica.GasMOP;

import etomica.space.Vector;
import etomica.space3d.Vector3D;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

public class XYZWriter {

    public static void writeXYZ(
            String filename,
            Map<Integer, String> atomMap,          // atomId -> "C_1" or "O_2" etc
            Map<Integer, Vector> positionMap,    // atomId -> position
            String comment
    ) throws IOException {

        // Sort by atomId so output is deterministic
        Map<Integer, String> sortedAtomMap = new TreeMap<>(atomMap);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filename))) {

            // Line 1: number of atoms
            writer.write(Integer.toString(sortedAtomMap.size()));
            writer.newLine();

            // Line 2: comment
            writer.write(comment == null ? "" : comment);
            writer.newLine();

            // Atom lines
            for (Map.Entry<Integer, String> e : sortedAtomMap.entrySet()) {
                int atomId = e.getKey();
                String typeName = e.getValue();          // e.g. "O_2" or "C_5"
                Vector pos = positionMap.get(atomId);

                if (pos == null) {
                    throw new IllegalArgumentException("Missing position for atomId=" + atomId);
                }

                // XYZ wants element symbol, not your numbered type. Convert "O_2" -> "O"
                String element = elementFromTypeName(typeName);

                String line = String.format("%-2s %12.6f %12.6f %12.6f",
                        element, pos.getX(0), pos.getX(1), pos.getX(2));

                writer.write(line);
                writer.newLine();
            }
        }
    }

    // "O_2" -> "O", "C_10" -> "C", "Cl_1" -> "Cl"
    private static String elementFromTypeName(String typeName) {
        if (typeName == null || typeName.isEmpty()) return "X";
        int underscore = typeName.indexOf('_');
        return (underscore > 0) ? typeName.substring(0, underscore) : typeName;
    }
}
