package etomica.virial.simulations.KnottedPolymer;

import etomica.config.IConformation;
import etomica.space.Space;
import etomica.space3d.Vector3D;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

/**
 * Conformation for star-polymers by reading the coordinates file.
 */
public class ConformationStarPolymerAll extends ConformationKnottedPolymerAll {

    protected double maxDiameter = 0.0;

    ConformationStarPolymerAll(Space space, String fileName, int index) {
        super(space, fileName, index);
    }

    @Override
    public HashMap<Integer, List<Vector3D>> readingFiles(String filePath) {
        int numOfConformation = 10000;
        HashMap<Integer, List<Vector3D>> starPolymers = new HashMap<>();
        FileInputStream inputStream = null;
        Scanner sc = null;
        try {
            inputStream = new FileInputStream(filePath);
            sc = new Scanner(inputStream, "UTF-8");

            if (sc.hasNextLine()) {
                for (int num = 0; num < numOfConformation; num++) {
                    String line = sc.nextLine();
                    List<Vector3D> coordinates = new ArrayList<>();
                    int nBead = Integer.parseInt(line);
                    sc.nextLine();  // Skip the second line (Atoms, TimeSteps)
                    for (int iBead = 0; iBead < nBead; iBead++) {
                        String[] coordStr = sc.nextLine().split("[ ]+");
                        double x = Double.parseDouble(coordStr[1]);
                        double y = Double.parseDouble(coordStr[2]);
                        double z = Double.parseDouble(coordStr[3]);
                        Vector3D coord = new Vector3D();
                        coord.E(x, y, z);
                        coordinates.add(coord);
                    }
                    starPolymers.put(num, coordinates);

                    double coreX = coordinates.get(nBead - 1).getX(0);
                    double coreY = coordinates.get(nBead - 1).getX(1);
                    double coreZ = coordinates.get(nBead - 1).getX(2);

                    double maxRadius = 0.0;
                    for (int index = 0; index < nBead - 1; index++) {
                        double beadX = coordinates.get(index).getX(0);
                        double beadY = coordinates.get(index).getX(1);
                        double beadZ = coordinates.get(index).getX(2);

                        double r2 = Math.pow(coreX - beadX, 2) +
                                Math.pow(coreY - beadY, 2) + Math.pow(coreZ - beadZ, 2);
                        if (r2 > maxRadius * maxRadius) {
                            maxRadius = Math.sqrt(r2);
                        }

                        maxDiameter = 2 * maxRadius;
                    }
                }

                if (sc.ioException() != null) {
                    throw sc.ioException();
                }
            }
        } catch (IOException e) {
            System.out.println(e.toString());
            System.out.println("Could not find file " + filePath);
        } finally {
            try {
                if (inputStream != null) {
                    inputStream.close();
                }
                if (sc != null) {
                    sc.close();
                }
            } catch (IOException e) {
                System.out.println(e.toString());
                System.out.println("Could not close file ");
            }
        }
        return starPolymers;
    }

    @Override
    public IConformation getNextConformation() {
        return null;
    }

    public double getMaxDiameter() {
        return maxDiameter;
    }
}
