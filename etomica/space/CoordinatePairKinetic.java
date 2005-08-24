package etomica.space;

import etomica.atom.AtomPair;

/*
 * History Created on Jan 24, 2005 by kofke
 */
public class CoordinatePairKinetic extends CoordinatePair {

    private final Vector dv;

    public CoordinatePairKinetic(Space space) {
        super(space);
        dv = space.makeVector();
    }
    
    public Vector resetV(AtomPair pair) {
        return resetV(pair.atom0.coord, pair.atom1.coord);
    }

    public Vector resetV(ICoordinate coord1, ICoordinate coord2) {
        c1 = (Coordinate)coord1;
        c2 = (Coordinate)coord2;
        return resetV();
    }

    public Vector resetV() {
        dv.Ev1Mv2(((ICoordinateKinetic) c2).velocity(),
                ((ICoordinateKinetic) c1).velocity());
        return dv;
    }

    public final Vector dv() {
        return dv;
    }

    public final double v2() {
        return dv.squared();
    }

    public final double vDotr() {
        return dr.dot(dv);
    }

    public void push(double impulse) {
        ((ICoordinateKinetic) c1).velocity().PEa1Tv1(+impulse, dr);
        ((ICoordinateKinetic) c2).velocity().PEa1Tv1(-impulse, dr);
    }
}