package etomica.parser.gromacs;

import etomica.space3d.Vector3D;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.StringJoiner;

public class GromacsGroParser {

    public static GroFile parseInput(Readable contents) {
        Scanner scanner = new Scanner(contents);
        String firstLine = scanner.nextLine(); // get time
        int numAtoms = scanner.nextInt();
        scanner.nextLine();
        List<GroLine> lines = new ArrayList<>();
        for (int i = 0; i < numAtoms; i++) {
            String line = scanner.nextLine();
            int resNum = Integer.parseInt(line.substring(0, 5).trim());
            String resName = line.substring(5, 10).trim();
            String atomName = line.substring(10, 15).trim();
            int atomNum = Integer.parseInt(line.substring(15, 20).trim());
            double px = Double.parseDouble(line.substring(20, 28).trim());
            double py = Double.parseDouble(line.substring(28, 36).trim());
            double pz = Double.parseDouble(line.substring(36, 44).trim());

            Vector3D vel = null;
            if (line.length() > 44) {
                double vx = Double.parseDouble(line.substring(44, 52).trim());
                double vy = Double.parseDouble(line.substring(52, 60).trim());
                double vz = Double.parseDouble(line.substring(60, 68).trim());
                vel = new Vector3D(vx, vy, vz);

                if (line.length() > 68) {
                    throw new RuntimeException("Malformed GRO line: " + line);
                }
            }
            lines.add(new GroLine(resNum, resName, atomName, atomNum, new Vector3D(px, py, pz), vel));
        }
        double boxX = scanner.nextDouble();
        double boxY = scanner.nextDouble();
        double boxZ = scanner.nextDouble();

        return new GroFile(lines, boxX, boxY, boxZ);
    }

    public static class GroFile {
        public final List<GroLine> lines;
        public final double boxX;
        public final double boxY;
        public final double boxZ;

        @Override
        public String toString() {
            return new StringJoiner(", ", GroFile.class.getSimpleName() + "[", "]")
                    .add("boxX=" + boxX)
                    .add("boxY=" + boxY)
                    .add("boxZ=" + boxZ)
                    .add("lines=" + lines)
                    .toString();
        }

        public GroFile(List<GroLine> lines, double boxX, double boxY, double boxZ) {
            this.lines = lines;
            this.boxX = boxX;
            this.boxY = boxY;
            this.boxZ = boxZ;
        }
    }

    public static class GroLine {
        public final int residueNumber;
        public final String residueName;
        public final String atomName;
        public final int atomNumber;
        public final Vector3D position;
        public final Vector3D velocity;

        public GroLine(int residueNumber, String residueName, String atomName, int atomNumber, Vector3D position, Vector3D velocity) {
            this.residueNumber = residueNumber;
            this.residueName = residueName;
            this.atomName = atomName;
            this.atomNumber = atomNumber;
            this.position = position;
            this.velocity = velocity;
        }

        @Override
        public String toString() {
            return new StringJoiner(", ", GroLine.class.getSimpleName() + "[", "]")
                    .add("residueNumber=" + residueNumber)
                    .add("residueName='" + residueName + "'")
                    .add("atomName='" + atomName + "'")
                    .add("atomNumber=" + atomNumber)
                    .add("position=" + position)
                    .add("velocity=" + velocity)
                    .toString();
        }
    }

    public static void main(String[] args) throws IOException {
        System.out.println(parseInput(Files.newBufferedReader(Paths.get("/home/alex/workspace/mosdef-test/ethane-box.gro"))));
    }
}
