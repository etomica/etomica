package etomica.space;

import etomica.Space;

/*
 * History Created on Jan 24, 2005 by kofke
 */
public class CoordinatePairKinetic extends CoordinatePair {

    private final Vector dv;

    public CoordinatePairKinetic(Space space) {
        super(space);
        dv = space.makeVector();
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