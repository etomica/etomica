package simulate;

import java.awt.*;
import java.beans.Beans;

public class PotentialSquareWell extends Potential {

  private double coreDiameter, coreDiameterSquared;
  private double wellDiameter, wellDiameterSquared;
  private double lambda; //wellDiameter = coreDiameter * lambda
  private double epsilon;
  private transient final double[] r12 = new double[Space.D];
  private transient final double[] v12 = new double[Space.D];
//  private transient PairInteraction pair = new PairInteraction();

  public PotentialSquareWell(double coreDiameter, double lambda, double epsilon) {
    setCoreDiameter(coreDiameter);
    setLambda(lambda);
    setEpsilon(epsilon);
  }

 /* public PairInteraction computePairInteraction(Molecule molecule1, Molecule molecule2) {
    MoleculeAtomic disk1 = (MoleculeAtomic)molecule1;
    MoleculeAtomic disk2 = (MoleculeAtomic)molecule2;        
    space.uEr1Mr2(pair.rij,disk2.r,disk1.r);
    pair.rSquared = Space.v1Dv2(pair.rij, pair.rij);
    pair.energy = (pair.rSquared < wellDiameterSquared) ? -epsilon : 0.0;
    return pair;
  }
 */
 
  public void bump(Atom atom1, Atom atom2) {
    double eps = 1.0e-6;
    space.uEr1Mr2(r12,atom2.r,atom1.r);  //use instance method   //r2-r1
    Space.uEa1Tv1Ma2Tv2(v12,atom2.rm,atom2.p,atom1.rm,atom1.p);  //v2-v1
    double r2 = Space.v1S(r12);
    double bij = Space.v1Dv2(v12, r12);
    // ke is kinetic energy due to components of velocity
    double reduced_m = 1.0/((1.0/atom1.rm+ 1.0/atom2.rm)*atom1.rm*atom2.rm);
    double ke = bij*bij*reduced_m/(2.0*r2);
    double s, r2New;
    if(2*r2 < (coreDiameterSquared+wellDiameterSquared)) {   // Hard-core collision
      s = 2.0*reduced_m*bij/r2;
      r2New = r2;
    }
    else {    // Well collision
      if(bij > 0.0) {         // Separating
	    if(ke < epsilon) {     // Not enough kinetic energy to escape
	       s = 2.0*reduced_m*bij/r2;
	       r2New = (1-eps)*wellDiameterSquared; 
	    }
	    else {                 // Escape
//	  s = (0.5*bij/r - Math.sqrt(0.5*(ke - epsilon)))/r;
	       s = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsilon/reduced_m))/r2;
	       r2New = (1+eps)*wellDiameterSquared;
	    }
      }
      else {                  // Approaching
//	s = (0.5*bij/r + Math.sqrt(0.5*(ke + epsilon)))/r;
	     s = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*epsilon/reduced_m))/r2;
	     r2New = (1-eps)*wellDiameterSquared; 
      }
    }

    Space.uPEa1Tv1(atom1.p, s, r12);
    Space.uMEa1Tv1(atom2.p, s, r12);
    if(r2New == r2) return;
    
    double ratio = atom1.rm/atom2.rm;  // (mass2/mass1)
    double delta = (Math.sqrt(r2New/r2) - 1.0)/(1+ratio);
    Space.uPEa1Tv1(atom1.r, -ratio*delta, r12);
    Space.uPEa1Tv1(atom2.r, delta, r12);
  }

//----------------------------------------------------------------------

  public double collisionTime(Atom atom1, Atom atom2) {
    double discr = 0.0;

    space.uEr1Mr2(r12,atom2.r,atom1.r);  //use instance method   //r2-r1
    Space.uEa1Tv1Ma2Tv2(v12,atom2.rm,atom2.p,atom1.rm,atom1.p);  //v2-v1
    double bij = Space.v1Dv2(r12,v12);                           //r12 . v12
    double r2 = Space.v1Dv2(r12,r12);
    double v2 = Space.v1Dv2(v12,v12);

    double tij = Double.MAX_VALUE;

    if(r2 < wellDiameterSquared) {  // Already inside wells

      if(r2 < coreDiameterSquared) {   // Inside core; collision now if approaching, at well if separating
        return (bij < 0) ? 0.0 : (-bij + Math.sqrt(bij*bij - v2 * ( r2 - wellDiameterSquared )))/v2;
      }

      if(bij < 0.0) {    // Check for hard-core collision
	    discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
	    if(discr > 0) {  // Hard cores collide next
	      tij = (-bij - Math.sqrt(discr))/v2;
	    }
	    else {           // Moving toward each other, but wells collide next
	      discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
	      tij = (-bij + Math.sqrt(discr))/v2;
	    }
      }
      else {           // Moving away from each other, wells collide next
	      discr = bij*bij - v2 * ( r2 - wellDiameterSquared );  // This is always > 0
	      tij = (-bij + Math.sqrt(discr))/v2;
      }
    }
    else {              // Outside wells; look for collision at well
      if(bij < 0.0) {
	    discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
	    if(discr > 0) {
	      tij = (-bij - Math.sqrt(discr))/v2;
	    }
      }
    }
    return tij;
  }
  
    public double energy(Atom atom1, Atom atom2) {
        double r2 = space.r1Mr2_S(atom1.r, atom2.r);
        return ( r2 < wellDiameterSquared) ? 
                    ((r2 < coreDiameterSquared) ? Double.MAX_VALUE : -epsilon) : 0.0;
    }
 

    public double getCoreDiameter() {return coreDiameter;}
    public void setCoreDiameter(double c) {
        coreDiameter = c;
        coreDiameterSquared = c*c;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }

    public double getLambda() {return lambda;}
    public void setLambda(double lam) {
        lambda = lam;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }
    
    public double getEpsilon() {return epsilon;}
    public void setEpsilon(double eps) {
        epsilon = eps;
    }

}
  