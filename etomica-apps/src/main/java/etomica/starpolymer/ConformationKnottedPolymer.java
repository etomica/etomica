package etomica.starpolymer;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;
import etomica.space3d.Vector3D;

import java.util.List;

/**
 * General Class for the conformation of knotted polymers by reading the '.xyz' file
 * May 7, 2018
 */
public abstract class ConformationKnottedPolymer implements IConformation {
    protected final Space space;
    public final String fileName;
    public List<Vector3D> coordinates;
    public double maxDiameter;
    public int nBead;

    ConformationKnottedPolymer(Space space, String fileName) {
        this.space = space;
        this.fileName = fileName;
//        readingFiles(this.fileName);
        this.coordinates = readingFiles(this.fileName);
        this.maxDiameter = getMaxDiameter();

    }


    /**
     * @param fileName
     * @return Coordinates of a polymer as a list of vector 3D
     */
    public abstract List<Vector3D> readingFiles(String fileName);

    private Vector3D getCoordinates(int indexBead) {
        return coordinates.get(indexBead);
    }

    // TODO: 5/7/18 How to get a random configuration from the file instead of only the first one
    public abstract IConformation getNextConformation();

    @Override
    public void initializePositions(IAtomList atomList) {
        boolean test = false;
        if (test) {
            int nBead = atomList.getAtoms().size();
            for (int iBead = 0; iBead < nBead; iBead++) {
                if (iBead == 0) {
                    atomList.get(0).getPosition().E(coordinates.get(200));
                } else {
                    IAtom bead = atomList.getAtoms().get(iBead);
                    bead.getPosition().E(getCoordinates(iBead - 1));
                }
            }
        } else {
            int nBead = atomList.getAtoms().size();
            for (int iBead = 0; iBead < nBead; iBead++) {
                IAtom bead = atomList.getAtoms().get(iBead);
                bead.getPosition().E(getCoordinates(iBead));
            }
        }
    }

    public double getMaxDiameter() {
        List<Vector3D> coordinates = this.coordinates;

        double coreX = coordinates.get(nBead - 1).getX(0);
        double coreY = coordinates.get(nBead - 1).getX(1);
        double coreZ = coordinates.get(nBead - 1).getX(2);

        for (int index = 0; index < nBead - 1; index++) {
            double beadX = coordinates.get(index).getX(0);
            double beadY = coordinates.get(index).getX(1);
            double beadZ = coordinates.get(index).getX(2);

            double r = Math.sqrt(Math.pow(coreX - beadX, 2) +
                    Math.pow(coreY - beadY, 2) + Math.pow(coreZ - beadZ, 2));
            if (2 * r > maxDiameter) {
                maxDiameter = 2 * r;
            }
        }


        return maxDiameter;
    }
}
