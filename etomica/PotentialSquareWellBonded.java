package simulate;

/**
 * Identical to square-well potential, but modifies bonding state of atoms with collision,
 * when appropriate.
 */
public class PotentialSquareWellBonded extends PotentialSquareWell {
    
  public PotentialSquareWellBonded() {
    super();
  }
  public PotentialSquareWellBonded(Simulation sim) {
    super(sim);
  }
  public PotentialSquareWellBonded(double coreDiameter, double lambda, double epsilon) {
    super(coreDiameter, lambda, epsilon);
  }
  
  public PotentialSquareWellBonded(Simulation sim, double coreDiameter, double lambda, double epsilon) {
    super(sim, coreDiameter, lambda, epsilon);
  }

/**
 * Implements collision dynamics between two square-well atoms.
 * Includes all possibilities involving collision of hard cores, and collision of wells
 * both approaching and diverging
 */
  public void bump(AtomPair pair) {
    double eps = 1.0e-6;
    double r2 = pair.r2();
    double bij = pair.vDotr();
    // ke is kinetic energy due to components of velocity
    double reduced_m = 1.0/(pair.atom1().rm() + pair.atom2().rm());
    dr.E(pair.dr());
    double ke = bij*bij*reduced_m/(2.0*r2);
    double s, r2New;
    if(2*r2 < (coreDiameterSquared+wellDiameterSquared)) {   // Hard-core collision
      lastCollisionVirial = 2.0*reduced_m*bij;
      r2New = r2;
    }
    else {    // Well collision
      if(bij > 0.0) {         // Separating
	    if(ke < epsilon) {     // Not enough kinetic energy to escape
	       lastCollisionVirial = 2.0*reduced_m*bij;
	       r2New = (1-eps)*wellDiameterSquared; 
	    }
	    else {                 // Escape
//	  s = (0.5*bij/r - Math.sqrt(0.5*(ke - epsilon)))/r;
	       lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsilon/reduced_m));
	       r2New = (1+eps)*wellDiameterSquared;
	       pair.atom1().atomLink[0] = null;
	       pair.atom2().atomLink[0] = null;
	    }
      }
      else {                  // Approaching
//	s = (0.5*bij/r + Math.sqrt(0.5*(ke + epsilon)))/r;
	     lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*epsilon/reduced_m));
	     r2New = (1-eps)*wellDiameterSquared;
	     pair.atom1().atomLink[0] = new Atom.Linker(pair.atom2());
	     pair.atom2().atomLink[0] = new Atom.Linker(pair.atom1());
      }
    }
    lastCollisionVirialr2 = lastCollisionVirial/r2;
    pair.cPair.push(lastCollisionVirialr2);
    if(r2New != r2) pair.cPair.setSeparation(r2New);
  }//end of bump
}//end of class
    