package simulate.space2D;

public final class CoordinatePair extends simulate.Space.CoordinatePair {  
    Coordinate c1;
    Coordinate c2;
    final Boundary boundary;
    final Vector dimensions;   //assumes this is not transferred between phases
    private final Vector dr = new Vector();  //note that dr is not cloned if this is cloned -- should change this if using dr in clones; also this makes cloned coordinatePairs not thread-safe
    private double drx, dry, dvx, dvy;
    public CoordinatePair() {boundary = new BoundaryNone(); dimensions = (Vector)boundary.dimensions();}
    public CoordinatePair(Space.Boundary b) {boundary = (Boundary)b; dimensions = (Vector)boundary.dimensions();}
    public void reset(Space.Coordinate coord1, Space.Coordinate coord2) {  //don't usually use this; instead set c1 and c2 directly, without a cast
        c1 = (Coordinate)coord1;
        c2 = (Coordinate)coord2;
        reset();
    }
    public void reset() {
        dr.x = c2.r.x - c1.r.x;
        dr.y = c2.r.y - c1.r.y;
        boundary.nearestImage(dr);
        drx = dr.x; 
        dry = dr.y;
        r2 = drx*drx + dry*dry;
        double rm1 = c1.parent().rm();
        double rm2 = c2.parent().rm();
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
        r2 = drx*drx + dry*dry;
    }
    public simulate.Space.Vector dr() {return dr;}
    public double dr(int i) {return (i==0) ? drx : dry;}
    public double dv(int i) {return (i==0) ? dvx : dvy;}
    public double v2() {
        return dvx*dvx + dvy*dvy;
    }
    public double vDotr() {
        return drx*dvx + dry*dvy;
    }
    public void push(double impulse) {  //changes momentum in the direction joining the atoms
        c1.p.x += impulse*drx;
        c1.p.y += impulse*dry;
        c2.p.x -= impulse*drx;
        c2.p.y -= impulse*dry;
    }
    public void setSeparation(double r2New) {
        double ratio = c2.parent().mass()*c1.parent().rm();  // (mass2/mass1)
        double delta = (Math.sqrt(r2New/this.r2()) - 1.0)/(1+ratio);
        c1.r.x -= ratio*delta*drx;
        c1.r.y -= ratio*delta*dry;
        c2.r.x += delta*drx;
        c2.r.y += delta*dry;
        //need call reset?
    }
}

