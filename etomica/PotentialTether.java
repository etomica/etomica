package simulate;

public class PotentialTether extends Potential implements PotentialHard {

  private double tetherLength, tetherLengthSquared;

  public PotentialTether() {
    setTetherLength(0.1);
  }

  public final double getTetherLength() {return tetherLength;}
  public final void setTetherLength(double t) {
      tetherLength = t;
      tetherLengthSquared = t*t;
  }

//----------------------------------------------------------------------

  public final void bump(AtomPair pair) {
        double r2 = pair.r2();
        double factor = 2.0/(pair.atom1().rm() + pair.atom2().rm())*pair.vDotr()/r2;
        pair.cPair.push(factor);
/*    parentPhase.space.uEr1Mr2(r12,atom2.r,atom1.r);  //use instance method   //r2-r1
    Space.uEa1Tv1Ma2Tv2(v12,atom2.rm,atom2.p,atom1.rm,atom1.p);  //v2-v1
    double r2 = tetherLengthSquared;
    double bij = Space.v1Dv2(v12, r12);
    double reduced_m = 1.0/((1.0/atom1.rm+ 1.0/atom2.rm)*atom1.rm*atom2.rm);
    double s = 2.0*reduced_m*bij/r2;  //same even if an inner-shell collision
    Space.uPEa1Tv1(atom1.p, s, r12);
    Space.uMEa1Tv1(atom2.p, s, r12); */
  }

//----------------------------------------------------------------------

  public final double collisionTime(AtomPair pair) {
    boolean minus;  //flag for +/- root
    double r2 = pair.r2();
    double bij = pair.vDotr();
 //       if(r2 < sig2) {return (bij > 0) ? Double.MAX_VALUE : 0.0;}  //inside wall; no collision
    if(r2 > tetherLengthSquared && bij >= 0) {return 0.0;}  //outside tether, moving apart; collide now
    double v2 = pair.v2();
    double discr = bij*bij - v2 * ( r2 - tetherLengthSquared );
    return (-bij + Math.sqrt(discr))/v2;
    }
  
    public double energy(AtomPair pair) {
        return (pair.r2() > tetherLengthSquared) ? Double.MAX_VALUE : 0.0;
    }
      
}
  