package etomica.virial;

import etomica.Space;
import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomPositionDefinitionSimple;
import etomica.exception.MethodNotImplementedException;
import etomica.space.ICoordinate;
import etomica.space.Vector;

public class CoordinatePairMolecular extends etomica.space.CoordinatePair {

    public CoordinatePairMolecular(Space space) {
        this(space,new AtomPositionDefinitionSimple());
    }
    
    public CoordinatePairMolecular(Space space, AtomPositionDefinition positionDefinition) {
        super(space);
        setAtomPositionDefinition(positionDefinition);
    }

    public Vector reset(AtomPair pair) {
        return reset(pair.atom0, pair.atom1);
    }

    public Vector reset(ICoordinate coord1, ICoordinate coord2) {
        throw new MethodNotImplementedException("can't reset a CoordinatePairMolecular with coordinates, silly!");
    }
    
    public Vector reset(Atom atom1, Atom atom2) {
        a1 = atom1;
        a2 = atom2;
        return reset();
    }

    public Vector reset() {
        dr.Ev1Mv2(apd.position(a2), apd.position(a1));
        nearestImageTransformer.nearestImage(dr);
        return dr;
    }

    public void nudge(double rDelta) {
        throw new MethodNotImplementedException("you really don't want to be here.  go away");
    }
    
    public void setAtomPositionDefinition(AtomPositionDefinition def) {
        apd = def;
    }
    
    protected Atom a1, a2;
    protected AtomPositionDefinition apd;
}