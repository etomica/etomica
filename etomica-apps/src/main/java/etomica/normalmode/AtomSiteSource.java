package etomica.normalmode;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculeAgentManager;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.OrientationFull3D;

public class AtomSiteSource implements AtomLeafAgentManager.AgentSource<Vector> {


    public AtomSiteSource(Space space) {
        this.space = space;
    }
    public Class getMoleculeAgentClass() {
        return Vector.class;
    }
    public Vector makeAgent(IAtom atom, Box box) {
        Vector agent = space.makeVector();
        agent.E(atom.getPosition());
        return agent;
    }
    public void releaseAgent(Vector agent, IAtom atom, Box box) {
        //nothing to do
    }

    private final Space space;
}