package etomica.space;

import etomica.NearestImageTransformer;
import etomica.Space;

/*
 * History Created on Jan 24, 2005 by kofke
 */
public class CoordinatePair {

    public CoordinatePair(Space space) {
        dr = space.makeVector();
    }

    public void setNearestImageTransformer(NearestImageTransformer b) {
        this.nearestImageTransformer = b;
    }

    public NearestImageTransformer getNearestImageTransformer() {
        return nearestImageTransformer;
    }

    public void reset(ICoordinate coord1, ICoordinate coord2) {
        c1 = (Coordinate)coord1;
        c2 = (Coordinate)coord2;
        reset();
    }

    public void reset() {
        dr.Ev1Mv2(c2.r, c1.r);
        nearestImageTransformer.nearestImage(dr);
    }

    public final double r2() {
        return dr.squared();
    }

    public final Vector dr() {
        return dr;
    }

    public void nudge(double rDelta) {
        c1.position().PEa1Tv1(-rDelta, dr);
        c2.position().PEa1Tv1(+rDelta, dr);
    }
    
    Coordinate c1;
    Coordinate c2;
    protected final Vector dr;
    protected NearestImageTransformer nearestImageTransformer = etomica.space.Boundary.NULL;

}