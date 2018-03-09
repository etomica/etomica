package etomica.normalmode;

import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculeAgentManager;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;

public class MoleculeSiteSourceNitrogen implements MoleculeAgentManager.MoleculeAgentSource {


    public MoleculeSiteSourceNitrogen(Space space, IMoleculePositionDefinition positionDefinition, MoleculeOrientationDefinition moleculeOrientationDefinition) {
        this.space = space;
        this.positionDefinition = positionDefinition;

        this.moleculeOrientationDefinition = moleculeOrientationDefinition;

    }

    public Class getMoleculeAgentClass() {
        return LatticeCoordinate.class;
    }

    public Object makeAgent(IMolecule molecule) {
        LatticeCoordinate latticeCoordinate = new LatticeCoordinate(space);
        latticeCoordinate.position.E(positionDefinition.position(molecule));
        latticeCoordinate.orientation.E(moleculeOrientationDefinition.getOrientation(molecule));
        return latticeCoordinate;
    }

    public void releaseAgent(Object agent, IMolecule molecule) {
        //nothing to do
    }

    private final Space space;
    protected final IMoleculePositionDefinition positionDefinition;
    protected final MoleculeOrientationDefinition moleculeOrientationDefinition;

    public static class LatticeCoordinate {
        public final Vector position;
        public final Orientation3D orientation;

        public LatticeCoordinate(Space space) {
            position = space.makeVector();
            orientation = new Orientation3D(space);
        }

    }

    public interface MoleculeOrientationDefinition {
        public IOrientation getOrientation(IMolecule molecule);

    }
}