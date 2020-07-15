package etomica.starpolymer;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;
import etomica.space3d.Vector3D;

import java.util.HashMap;
import java.util.List;

/**
 * General Class for the conformation of knotted polymers by reading the '.xyz' file
 * May 7, 2018
 */
public abstract class ConformationKnottedPolymerAll implements IConformation {
    protected final Space space;
    public final String fileName;
    public List<Vector3D> coordinates;
    public HashMap<Integer, List<Vector3D>> starPolymers;

    ConformationKnottedPolymerAll(Space space, String fileName, int index) {
        this.space = space;
        this.fileName = fileName;
        this.starPolymers = readingFiles(this.fileName);
        this.coordinates = starPolymers.get(index);
    }

    /**
     * @param fileName
     * @return Coordinates of a polymer as a list of vector 3D
     */
    public abstract HashMap<Integer, List<Vector3D>> readingFiles(String fileName);

    private Vector3D getCoordinates(int indexBead) {
        return coordinates.get(indexBead);
    }

    // TODO: 5/7/18 How to get a random configuration from the file instead of only the first one
    public abstract IConformation getNextConformation();

    @Override
    public void initializePositions(IAtomList atomList) {

        if (true) {
            int nBead = atomList.getAtoms().size();
            for (int iBead = 0; iBead < nBead; iBead++) {
                if (iBead == 0) {
                    IAtom bead = atomList.getAtoms().get(iBead);
                    bead.getPosition().E(getCoordinates(200));
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
}
