package etomica.virial;

import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomPositionDefinitionSimple;
import etomica.exception.MethodNotImplementedException;
import etomica.space.ICoordinate;
import etomica.space.Space;
import etomica.space.Vector;

public class CoordinatePairMolecular extends etomica.space.CoordinatePair {

    public CoordinatePairMolecular(Space space) {
        this(space,new AtomPositionDefinitionSimple());
    }
    
    public CoordinatePairMolecular(Space space, AtomPositionDefinition positionDefinition) {
        super(space);
        setAtomPositionDefinition(positionDefinition);
    }

    public void reset(AtomPair pair) {
        reset(pair.atom0, pair.atom1);
    }

    public void reset(ICoordinate coord1, ICoordinate coord2) {
        throw new MethodNotImplementedException("can't reset a CoordinatePairMolecular with coordinates, silly!");
    }
    
    public void reset(Atom atom1, Atom atom2) {
        a1 = atom1;
        a2 = atom2;
        reset();
    }

    public void reset() {
        dr.Ev1Mv2(apd.position(a2), apd.position(a1));
        nearestImageTransformer.nearestImage(dr);
    }

    public void setAtomPositionDefinition(AtomPositionDefinition def) {
        apd = def;
    }
    
    protected Atom a1, a2;
    protected AtomPositionDefinition apd;
}