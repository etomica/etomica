package etomica.space1d;

import etomica.NearestImageTransformer;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public final class CoordinatePair extends etomica.space.CoordinatePair {
    Coordinate c1;
    Coordinate c2;
    private final Vector dr = new Vector(); //note that dr is not cloned if this is cloned -- this should be changed if cloned vectors use dr; also this makes cloned coordinatePairs not thread-safe
    private double drx, dvx;
	private NearestImageTransformer nearestImageTransformer = etomica.space.Boundary.NULL;

	public void setNearestImageTransformer(NearestImageTransformer b) {this.nearestImageTransformer = b;}
	public NearestImageTransformer getNearestImageTransformer() {return nearestImageTransformer;}		
    public double r2() {return drx*drx;}
    public void reset(etomica.space.Coordinate coord1, etomica.space.Coordinate coord2) {  //don't usually use this; instead set c1 and c2 directly, without a cast
        c1 = (Coordinate)coord1;
        c2 = (Coordinate)coord2;
        reset();
    }
    public void trueReset(etomica.space.Coordinate coord1, etomica.space.Coordinate coord2, double falseTime) {
        c1 = (Coordinate)coord1;
        c2 = (Coordinate)coord2;
        trueReset(falseTime);
    }
    public void reset() {
        dr.x = c2.r.x - c1.r.x;
        c1.atom.node.parentPhase().boundary().nearestImage(dr);
        drx = dr.x;
    }
    public void trueReset(double falseTime) {
        resetV();
        dr.Ev1Mv2(c2.r,c1.r);

        dr.x += falseTime * dvx;
        nearestImageTransformer.nearestImage(dr);
        drx = dr.x;
    }
    public void resetV() {
        double rm1 = c1.rm();
        double rm2 = c2.rm();
        dvx = (rm2*c2.p.x - rm1*c1.p.x);  
    }
    /**
     * Recomputes pair separation, with atom 2 shifted by the given vector
     * Does not apply any PBC, regardless of boundary chosen for space
     */
    public void reset(Vector M) {
        dr.x = c2.r.x - c1.r.x + M.x;
        drx = dr.x;
    }
    public etomica.space.Vector dr() {return dr;}
    public double dr(int i) {return drx;}
    public double dv(int i) {return dvx;}
    public double v2() {
        return dvx*dvx;
    }
    public double vDot(etomica.space.Vector u) {return vDot((Vector)u);}
    public double vDot(Vector u) {return dvx*u.x;}
    public double vDotr() {
        return drx*dvx;
    }
    public void push(double impulse) {  //changes momentum in the direction joining the atoms
        c1.p.x += impulse*drx;
        c2.p.x -= impulse*drx;
    }
    public void truePush(etomica.space.Vector u, double falseTime) {
        c1.p.PE(u);
        c2.p.ME(u);
        
        c1.r.PEa1Tv1(-falseTime*c1.rm(),u);
        c2.r.PEa1Tv1(falseTime*c2.rm(),u);
    }
    public void nudge(double rDelta) {
        double ratio = c2.mass()*c1.rm()*rDelta;
        c1.r.x -= ratio*dr.x;
        c2.r.x += ratio*dr.x;
    }
    public void setSeparation(double r2New) {
        double ratio = c2.mass()*c1.rm();  // (mass2/mass1)
        double delta = (Math.sqrt(r2New/this.r2()) - 1.0)/(1+ratio);
        c1.r.x -= ratio*delta*drx;
        c2.r.x += delta*drx;
        //need call reset?
    }
}
