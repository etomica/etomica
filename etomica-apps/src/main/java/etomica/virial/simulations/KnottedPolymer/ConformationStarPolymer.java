package etomica.virial.simulations.KnottedPolymer;

import etomica.config.IConformation;
import etomica.space.Space;
import etomica.space3d.Vector3D;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Conformation for star-polymers by reading the coordinates file.
 */
public class ConformationStarPolymer extends ConformationKnottedPolymer {

    ConformationStarPolymer(Space space, String fileName) {
        super(space, fileName);
    }

    @Override
    public List<Vector3D> readingFiles(String fileName) {
        List<Vector3D> coordinates = new ArrayList<>();
        FileReader fileReader;
        try {
            fileReader = new FileReader(fileName);
        } catch (IOException e) {
            throw new RuntimeException("Cannot open " + fileName + ", caught IOException: " + e.getMessage());
        }

        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            this.nBead = Integer.parseInt(bufReader.readLine());
            int nBead = this.nBead;

            // skip the first and read the second instead
//            int ok = 203;
//            while (ok > 0){
//                bufReader.readLine();
//                ok--;
//            }
            bufReader.readLine();  // Skip the second line (Atoms, TimeSteps)
            for (int iBead = 0; iBead < nBead; iBead++) {
                String[] coordStr = bufReader.readLine().split("[ ]+");
                double x = Double.parseDouble(coordStr[1]);
                double y = Double.parseDouble(coordStr[2]);
                double z = Double.parseDouble(coordStr[3]);
                Vector3D coord = new Vector3D();
                coord.E(x, y, z);
                coordinates.add(coord);
            }

            fileReader.close();
        } catch (IOException e1) {
            e1.printStackTrace();
        }

        return coordinates;
    }

    @Override
    public IConformation getNextConformation() {
        return null;
    }

}
