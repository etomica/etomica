package etomica.modules.reactionequilibrium;
import etomica.Atom;
import etomica.AtomPair;
import etomica.Space;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialReactive;
import etomica.space.CoordinatePairKinetic;
import etomica.units.Dimension;

/**
 * Similar to square-well potential, but considers and alters bonding states with collisions.
 * Each atom may bind, in the form of the square-well attraction, to one other atom.  Two atoms
 * approaching each other with this potential may interact as follows:
 * <ul>
 *   <li>If neither is bound, or if they are bound to each other, 
 *       they will interact as square-well atoms.
 *   <li>If   one or both is bound, each to another atom, they will act as hard
 * spheres of diameter equal to the well diameter.
 * </ul>
 * The potential is similar to P2SquareWellBondedBarrier, but there is no
 * accounting for a barrier, and no possibility for one atom to dislodge the
 * bonding partner of another directly in a single collision.
 *
 * @author David Kofke
 */
public class P2SquareWellBonded extends P2SquareWell implements PotentialReactive, Atom.AgentSource {
  
  public final int idx;
  
  private PotentialReactive.BondChangeData[] bondData;
  private double barrier;
  private final CoordinatePairKinetic cPair;
  
  public P2SquareWellBonded(Space parent,int idx, double coreDiameter, double lambda, double epsilon) {
	super(parent, coreDiameter, lambda, epsilon);
	cPair = (CoordinatePairKinetic) parent.makeCoordinatePair();
  //constructs 2 BondChangeData classes, one for each atom in possible bondchange event
	this.idx = Atom.requestAgentIndex(this);
	bondData = new PotentialReactive.BondChangeData[2];
	for(int i=0; i<2; i++) {
		bondData[i] = new PotentialReactive.BondChangeData();
		bondData[i].oldPartners = new Atom[1];
		bondData[i].newPartners = new Atom[1];
	}
  }
 /* public P2SquareWellBonded(int idx, double coreDiameter, double lambda, double epsilon) {
  	this(Simulation.instance.hamiltonian.potential, idx, coreDiameter, lambda, epsilon);
  }
  public P2SquareWellBonded(PotentialGroup parent, int idx, double coreDiameter, double lambda, double epsilon) {
    super(parent, coreDiameter, lambda, epsilon);
  //constructs 2 BondChangeData classes, one for each atom in possible bondchange event
    this.idx = idx;
    bondData = new PotentialReactive.BondChangeData[2];
    for(int i=0; i<2; i++) {
        bondData[i] = new PotentialReactive.BondChangeData();
        bondData[i].oldPartners = new Atom[1];
        bondData[i].newPartners = new Atom[1];
    }
  }
  */
  /**
   * Implementation of Atom.AgentSource interface, returning null.  Agent
   * in atom is used to hold bonding partner.
   * @param a ignored
   * @return Object always null
   */
  public Object makeAgent(Atom a) {
	return null;
  }
   

  public void setBarrier(double b) {barrier = b;}
  public double getBarrier() {return barrier;}
  public Dimension getBarrierDimension() {return Dimension.ENERGY;}
  
  public PotentialReactive.BondChangeData[] getBondChangeData() {return bondData;}

/**
 * Computes next time of collision of the two atoms, assuming free-flight kinematics.
 */
  public double collisionTime(AtomPair pair, double falsetime) {

    Atom a1Partner = (Atom)pair.atom0.allatomAgents[idx];
    cPair.reset(pair);
    //inside well but not mutually bonded;  collide now if approaching
    if(a1Partner != pair.atom1 && cPair.r2() < wellDiameterSquared) 
        return (cPair.vDotr() < 0.0) ? 0.0 : Double.MAX_VALUE;
    //mutually bonded, or outside well; collide as SW
    else return super.collisionTime(pair, falsetime);
  }

/**
 * Implements collision dynamics between the atom pair.
 */
  public void bump(AtomPair pair) {
    cPair.reset(pair);
  	double nudge = 0;
    double eps = 1.0e-6;
    double r2 = cPair.r2();
    double bij = cPair.vDotr();
    // ke is kinetic energy due to components of velocity
    double reduced_m = 1.0/(pair.atom0.type.rm() + pair.atom1.type.rm());
    dr.E(cPair.dr());
    double ke = bij*bij*reduced_m/(2.0*r2);
    double s, r2New;
    
    Atom a1 = pair.atom0;
    Atom a2 = pair.atom1;
    Atom a1Partner = (Atom)a1.allatomAgents[idx];
    Atom a2Partner = (Atom)a2.allatomAgents[idx];
    
    if(a1Partner == a2) //atoms are bonded to each other
        if(2*r2 < (coreDiameterSquared+wellDiameterSquared)) {   // Hard-core collision
            lastCollisionVirial = 2.0*reduced_m*bij;
            r2New = r2;
            nudge = eps;
            bondData[0].atom = null;
        }
        else {    // Well collision; assume separating because mutually bonded
	        if(ke < epsilon) {     // Not enough kinetic energy to escape
	            lastCollisionVirial = 2.0*reduced_m*bij;
	            r2New = (1-eps)*wellDiameterSquared; 
	            bondData[0].atom = null;
	            nudge = -eps;
	        }
	        else {                 // Escape
	            lastCollisionVirial = reduced_m*bij - Math.sqrt(2.*reduced_m*r2*(ke-epsilon));
	            r2New = (1+eps)*wellDiameterSquared;
	            bondData[0].atom = a1;
	            bondData[0].oldPartners[0] = a2;
	            bondData[0].newPartners[0] = null;
	            bondData[1].atom = a2;
	            bondData[1].oldPartners[0] = a1;
	            bondData[1].newPartners[0] = null;
	   //         bondData[2].atom = null;
	            a1.allatomAgents[idx] = null;
	            a2.allatomAgents[idx] = null;
	            nudge = eps;
	        }//end if(ke < epsilon)
        }//end if(2*r2...
    else {                  //not bonded to each other
        if(r2 < (1-eps)*wellDiameterSquared) {   //already inside well; push away
            lastCollisionVirial = 2.0*reduced_m*bij;
            r2New = r2;
            bondData[0].atom = null;
            nudge = eps;
        }
        else { //well collision; decide whether to bond or have hard repulsion
        	if(a1Partner != null || a2Partner != null) { //at least one is already bound; repel
				lastCollisionVirial = 2.0*reduced_m*bij;
				r2New = (1-eps)*wellDiameterSquared;
				bondData[0].atom = null;
				nudge = eps;
        	} else { //neither is taken; bond to each other
				lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*epsilon/reduced_m));
				r2New = (1-eps)*wellDiameterSquared;
				a1.allatomAgents[idx] = a2;
				a2.allatomAgents[idx] = a1; 
				bondData[0].atom = a1;
			    bondData[0].newPartners[0] = a2;
			    bondData[0].oldPartners[0] = null;
			
			    bondData[1].atom = a2;
			    bondData[1].newPartners[0] = a1;
			    bondData[1].oldPartners[0] = null;
			    nudge = -eps;
			}
        

        }//end if(r2 < wellDiameterSquared) else
    } //end if(a1Partner == a2) else
    lastCollisionVirialr2 = lastCollisionVirial/r2;
    cPair.push(lastCollisionVirialr2);
    cPair.nudge(nudge);
  }//end of bump
}//end of class
    