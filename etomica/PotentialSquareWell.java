package simulate;

import java.awt.*;
import java.beans.Beans;

public class PotentialSquareWell extends Potential implements PotentialHard {

  private double coreDiameter, coreDiameterSquared;
  private double wellDiameter, wellDiameterSquared;
  private double lambda; //wellDiameter = coreDiameter * lambda
  private double epsilon;

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
 
  public void bump(AtomPair pair) {
    double eps = 1.0e-6;
    double r2 = pair.r2();
    double bij = pair.vDotr();
    // ke is kinetic energy due to components of velocity
    double reduced_m = 1.0/(pair.atom1().rm() + pair.atom2().rm());
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

    pair.cPair.push(s);
    if(r2New != r2) pair.cPair.setSeparation(r2New);
  }

//----------------------------------------------------------------------

  public double collisionTime(AtomPair pair) {
    double discr = 0.0;

    double bij = pair.vDotr();
    double r2 = pair.r2();
    double v2 = pair.v2();

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
  
    public double energy(AtomPair pair) {
        double r2 = pair.r2();
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
  