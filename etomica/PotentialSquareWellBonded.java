package etomica;

/**
 * Identical to square-well potential, but modifies bonding state of atoms with collision,
 * when appropriate.
 */
public class PotentialSquareWellBonded extends PotentialSquareWell implements Potential.Reactive {
  
  private Potential.Reactive.BondChangeData[] bondData;
  
  public PotentialSquareWellBonded() {
    super();
    makeBondData();
  }
  public PotentialSquareWellBonded(Simulation sim) {
    super(sim);
    makeBondData();
  }
  public PotentialSquareWellBonded(double coreDiameter, double lambda, double epsilon) {
    super(coreDiameter, lambda, epsilon);
    makeBondData();
  }
  
  public PotentialSquareWellBonded(Simulation sim, double coreDiameter, double lambda, double epsilon) {
    super(sim, coreDiameter, lambda, epsilon);
    makeBondData();
  }
  
  //constructs 3 BondChangeData classes, one for each atom in possible bondchange event
  private void makeBondData() {
    bondData = new Potential.Reactive.BondChangeData[3];
    for(int i=0; i<3; i++) {
        bondData[i] = new Potential.Reactive.BondChangeData();
        bondData[i].oldPartners = new Atom[1];
        bondData[i].newPartners = new Atom[1];
    }
  }
  
  public Potential.Reactive.BondChangeData[] getBondChangeData() {return bondData;}

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
      bondData[0].atom = null;
    }
    else {    // Well collision
      if(bij > 0.0) {         // Separating
	    if(ke < epsilon) {     // Not enough kinetic energy to escape
	       lastCollisionVirial = 2.0*reduced_m*bij;
	       r2New = (1-eps)*wellDiameterSquared; 
	       bondData[0].atom = null;
	    }
	    else {                 // Escape
//	  s = (0.5*bij/r - Math.sqrt(0.5*(ke - epsilon)))/r;
	       lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsilon/reduced_m));
	       r2New = (1+eps)*wellDiameterSquared;
	       bondData[0].atom = pair.atom1();
	       bondData[0].oldPartners[0] = pair.atom2();
	       bondData[0].newPartners[0] = null;
	       bondData[1].atom = pair.atom2();
	       bondData[1].oldPartners[0] = pair.atom1();
	       bondData[1].newPartners[0] = null;
	       bondData[2].atom = null;
	       pair.atom1().atomLink[0] = null;
	       pair.atom2().atomLink[0] = null;
	    }
      }
      else {                  // Approaching
//	s = (0.5*bij/r + Math.sqrt(0.5*(ke + epsilon)))/r;
        lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*epsilon/reduced_m));
        r2New = (1-eps)*wellDiameterSquared;
        
        bondData[0].atom = pair.atom1();
        bondData[0].newPartners[0] = pair.atom2();
        if(pair.atom1().atomLink[0] != null) 
            bondData[0].oldPartners[0] = pair.atom1().atomLink[0].atom();
        else
            bondData[0].oldPartners[0] = null;
        bondData[1].atom = pair.atom2();
        bondData[1].newPartners[0] = pair.atom1();
        if(pair.atom2().atomLink[0] != null) 
            bondData[1].oldPartners[0] = pair.atom2().atomLink[0].atom();
        else
            bondData[1].oldPartners[0] = null;
        bondData[2].atom = null;
        pair.atom1().atomLink[0] = new Atom.Linker(pair.atom2());
        pair.atom2().atomLink[0] = new Atom.Linker(pair.atom1());
      }
    }
    lastCollisionVirialr2 = lastCollisionVirial/r2;
    pair.cPair.push(lastCollisionVirialr2);
    if(r2New != r2) pair.cPair.setSeparation(r2New);
  }//end of bump
}//end of class
    