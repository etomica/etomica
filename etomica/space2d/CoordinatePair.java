package etomica.space2d;

import etomica.NearestImageTransformer;
import etomica.space.Coordinate;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public final class CoordinatePair extends CoordinatePair {
    Coordinate c1;
    Coordinate c2;
    private final Vector dr = new Vector();  //note that dr is not cloned if this is cloned -- should change this if using dr in clones; also this makes cloned coordinatePairs not thread-safe
    private double drx, dry, dvx, dvy;
	private NearestImageTransformer nearestImageTransformer = Boundary;
    public double r2() {return drx*drx + dry*dry;}
	public void setNearestImageTransformer(NearestImageTransformer b) {this.nearestImageTransformer = b;}
	public NearestImageTransformer getNearestImageTransformer() {return nearestImageTransformer;}		
    public void reset(Coordinate coord1, Coordinate coord2) {  //don't usually use this; instead set c1 and c2 directly, without a cast
        c1 = (Coordinate)coord1;
        c2 = (Coordinate)coord2;
        reset();
    }
    public void trueReset(Coordinate coord1, Coordinate coord2, double falseTime) {
        c1 = (Coordinate)coord1;
        c2 = (Coordinate)coord2;
        trueReset(falseTime);
    }
    public void reset() {
        dr.x = c2.r.x - c1.r.x;
        dr.y = c2.r.y - c1.r.y;
        nearestImageTransformer.nearestImage(dr);
        drx = dr.x; 
        dry = dr.y;
    }
    public void trueReset(double falseTime) {
        resetV();
        dr.Ev1Mv2(c2.r,c1.r);

        dr.x += falseTime * dvx;
        dr.y += falseTime * dvy;
        nearestImageTransformer.nearestImage(dr);
        drx = dr.x;
        dry = dr.y;
    }
    public void resetV() {
        double rm1 = c1.rm();
        double rm2 = c2.rm();
        dvx = (rm2*c2.p.x - rm1*c1.p.x);  
        dvy = (rm2*c2.p.y - rm1*c1.p.y);  
    }
    /**
     * Recomputes pair separation, with atom 2 shifted by the given vector
     * Does not apply any PBC, regardless of boundary chosen for space
     */
    public void reset(Vector M) {
        dr.x = c2.r.x - c1.r.x + M.x;
        dr.y = c2.r.y - c1.r.y + M.y;
        drx = dr.x;
        dry = dr.y;
    }
    public Vector dr() {return dr;}
    public double dr(int i) {return (i==0) ? drx : dry;}
    public double dv(int i) {return (i==0) ? dvx : dvy;}
    public double v2() {
        return dvx*dvx + dvy*dvy;
    }
    public double vDot(Vector u) {return vDot((Vector)u);}
    public double vDot(Space2D.Vector u) {return dvx*u.x + dvy*u.y;}
    public double vDotr() {
        return drx*dvx + dry*dvy;
    }
    public void push(double impulse) {  //changes momentum in the direction joining the atoms
        c1.p.x += impulse*drx;
        c1.p.y += impulse*dry;
        c2.p.x -= impulse*drx;
        c2.p.y -= impulse*dry;
    }
    public void truePush(Vector u, double falseTime) {
        c1.p.PE(u);
        c2.p.ME(u);
        
        c1.r.PEa1Tv1(-falseTime*c1.rm(),u);
        c2.r.PEa1Tv1(falseTime*c2.rm(),u);
     }
    public void nudge(double rDelta) {
        double ratio = c2.mass()*c1.rm()*rDelta;
        c1.r.x -= ratio*dr.x;
        c1.r.y -= ratio*dr.y;
        c2.r.x += ratio*dr.x;
        c2.r.y += ratio*dr.y;
    }
    public void setSeparation(double r2New) {
        double ratio = c2.mass()*c1.rm();  // (mass2/mass1)
        double delta = (Math.sqrt(r2New/this.r2()) - 1.0)/(1+ratio);
        c1.r.x -= ratio*delta*drx;
        c1.r.y -= ratio*delta*dry;
        c2.r.x += delta*drx;
        c2.r.y += delta*dry;
        //need call reset?
    }
}