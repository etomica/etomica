package etomica.space;
import etomica.Space;

public final class CoordinatePair extends Space.CoordinatePair {
    Coordinate c1;
    Coordinate c2;
    //put this back 
    private final Space.Vector dr;
    private final Space.Vector dv;
    //in place of this
 //   private final etomica.space.continuum.Vector3D dr;
 //   private final etomica.space.continuum.Vector dv;
    
    public CoordinatePair(Space space) {
        super();
        dr = space.makeVector();
        dv = space.makeVector();
  //      dr = (etomica.space.continuum.Vector3D)space.makeVector();
  //      dv = (etomica.space.continuum.Vector)space.makeVector();
    }

    public void reset(Space.Coordinate coord1, Space.Coordinate coord2) {
        c1 = (Coordinate)coord1;
        c2 = (Coordinate)coord2;
        reset();
    }
    public void reset() {
//        System.out.println("inside coordintepair reset");
        dr.Ev1Mv2(c2.r, c1.r);
        c1.atom.node.parentPhase().boundary().nearestImage(dr);
        r2 = dr.squared();
        double rm1 = c1.rm();
        double rm2 = c2.rm();
        dv.Ea1Tv1(c2.rm(), c2.p);
        dv.PEa1Tv1(-c1.rm(), c1.p); //dv = p2/m2 - p1/m1;
    }
            
    public void reset(Space.Vector M) {
        dr.Ev1Mv2(c2.r, c1.r);
        dr.PE(M);
        r2 = dr.squared();
    }
    public double r2() {
        return r2;
      //  dr.Ev1Mv2(c2.r, c1.r);
      //  c1.atom.node.parentPhase().boundary().nearestImage(dr);
      //  return dr.squared();
    }
            
    public Space.Vector dr() {return dr;}
    public double dr(int i) {return dr.component(i);}
    public double dv(int i) {return dv.component(i);}
    public double v2() {return dv.squared();}
    public double vDot(Space.Vector u) {return dv.dot(u);}
    public double vDotr() {return dv.dot(dr);}
    public void push(double impulse) {
        c1.p.PEa1Tv1(+impulse, dr);
        c2.p.PEa1Tv1(-impulse, dr);
        //c1.p.x += impulse*dr.x;
        //c2.p.x -= impulse*dr.x;
    }
    public void setSeparation(double r2New) {
        double ratio = c2.mass()*c1.rm();
        double delta = (Math.sqrt(r2New/this.r2()) - 1.0)/(1 + ratio);
        c1.r.PEa1Tv1(-ratio*delta, dr);
        c2.r.PEa1Tv1(+ratio*delta, dr);
        //c1.r.x -= ratio*delta*dr.x;
        //c2.r.x += ratio*delta*dr.x;
    }
}//end of CoordinatePair
        
