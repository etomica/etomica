package etomica.modules.reactionequilibrium;
import etomica.Default;
import etomica.Space;
import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.potential.P2SquareWell;
import etomica.space.CoordinatePairKinetic;
import etomica.space.ICoordinateKinetic;
import etomica.space.Vector;
import etomica.units.Dimension;


/**
 * Similar to square-well potential, but considers and alters bonding states
 * with collisions. Each atom may bind, in the form of the square-well
 * attraction, to one other atom. Two atoms approaching each other with this
 * potential may interact as follows:
 * <ul>
 * <li>If neither is bound, or if they are bound to each other, they will
 * interact as square-well atoms.
 * <li>If one or both is bound, each to another atom, they will act as hard
 * spheres of diameter equal to the well diameter.
 * </ul>
 * The potential is similar to P2SquareWellBondedBarrier, but there is no
 * accounting for a barrier, and no possibility for one atom to dislodge the
 * bonding partner of another directly in a single collision.
 * 
 * @author David Kofke
 */
public class P2SquareWellBonded extends P2SquareWell {

	// This is an  index for an array Array ==> {1,2,3,4,5,6...}
	//			          idx = 2      | 
	public final int idx;
	private double barrier;

	// *** Below are different types of constructor functions that Make P2SqaredWells ***	

//	public P2SquareWellBonded(Space space, double coreDiameter, double lambda, double epsilon) {
//		super(space, coreDiameter, lambda, epsilon);
		//constructs 2 BondChangeData classes, one for each atom in possible
		// bondchange event
//	}

	public P2SquareWellBonded(Space space, int idx, double coreDiameter,double lambda, double epsilon) {
		super(space, coreDiameter, lambda, epsilon);
		this.idx = idx;
	}

	/**
	 * Implementation of Atom.AgentSource interface, returning null. Agent in atom is used to hold bonding partner.
	 * @param a ignored  
	 * @return Object always null
	 */

	// **** This makeAgent function adds a peice of data to the atom class, we are adding 
	//	An array of 2 atoms, so from every atom in our simulation there is an array [1,2]
	// 	in each space is saved an atom which "points" to another atom in the simulation which
	//	That atom is bonded to. 
	
	/*
	public boolean full(Atom a) {

	
	}
	*/

	public Object makeAgent(Atom a) {
		return new Atom[2];
	}


	public void setBarrier(double b) {
		barrier = b;
	}

	public double getBarrier() {
		return barrier;
	}

	public Dimension getBarrierDimension() {
		return Dimension.ENERGY;
	}

	/**
	 * Computes next time of collision of the two atoms, assuming free-flight
	 * kinematics.
	 */

	public double collisionTime(AtomSet atoms, double falseTime) {
		
		if (Default.FIX_OVERLAP) {
			
			// ** Makes 2 things, and atomPair pair, 
			AtomPair pair = (AtomPair) atoms;
			Atom a1Partner = (Atom) pair.atom0.allatomAgents[idx];
			
			cPair.reset(pair);
			((CoordinatePairKinetic) cPair).resetV();
			dr.E(cPair.dr());
			Vector dv = ((CoordinatePairKinetic) cPair).dv();
			dr.PEa1Tv1(falseTime, dv);
			double r2 = dr.squared();
			double bij = dr.dot(dv);

			//inside well but not mutually bonded; collide now if approaching
			if ((a1Partner != pair.atom1 && r2 < wellDiameterSquared)
					|| (a1Partner == pair.atom1 && r2 < coreDiameterSquared))
				return (bij < 0.0) ? falseTime : Double.POSITIVE_INFINITY;
		}
		//mutually bonded, or outside well; collide as SW
		double time = super.collisionTime(atoms, falseTime);
//		if(!Double.isInfinite(time)) System.out.println("Collision time: "+time+" for "+atoms.toString());
		return time;
	}

	
	public void bump(AtomSet atoms, double falseTime) {
		
		// *** Data Declaration Section

		AtomPair pair = (AtomPair) atoms;
		cPair.reset(pair);
		((CoordinatePairKinetic) cPair).resetV();
		dr.E(cPair.dr());
		Vector dv = ((CoordinatePairKinetic) cPair).dv();
		dr.PEa1Tv1(falseTime, dv);
		double r2 = dr.squared();
		double bij = dr.dot(dv);
		double nudge = 0;
		double eps = 1.0e-10;
		
		// ke is kinetic energy due to components of velocity
		Atom a0 = pair.atom0;
		Atom a1 = pair.atom1;
        double rm0 = ((AtomTypeLeaf)a0.type).rm();
        double rm1 = ((AtomTypeLeaf)a1.type).rm();
//		System.out.println("Bumping "+pair.toString());
		double reduced_m = 2.0 / (rm0 + rm1);
		double ke = bij * bij * reduced_m / (4.0 * r2);
		
		Atom a0Partner1 = ((Atom[]) a0.allatomAgents[idx])[0];
		Atom a0Partner2 = ((Atom[]) a0.allatomAgents[idx])[1];
		Atom a1Partner1 = ((Atom[]) a1.allatomAgents[idx])[0];
		Atom a1Partner2 = ((Atom[]) a1.allatomAgents[idx])[1];

		boolean a0Saturated = (a0Partner1 != null) && (a0Partner2 != null);
		boolean a1Saturated = (a1Partner1 != null) && (a1Partner2 != null);

		if (a0Partner1 == a1 || a0Partner2 == a1) {	//atoms are bonded to each
            if (a1Partner1 != a0 && a1Partner2 != a0) {
                throw new IllegalStateException("partner mismatch");      
            }

			if (2 * r2 < (coreDiameterSquared + wellDiameterSquared)) { // Hard-core collision															
				lastCollisionVirial = reduced_m * bij;
				//System.out.println("core");
	
			} else { 				// Well collision assume separating because mutually bonded

				if (ke < epsilon) 		// Not enough kinetic energy to escape
				{ 
					lastCollisionVirial = reduced_m * bij;
					nudge = -eps;
					//System.out.println("no escape");
				} 
				else 				// Enough energy to Escape, atoms escape
				{ 				
					
				lastCollisionVirial = 0.5 * reduced_m * bij- Math.sqrt(reduced_m * r2 * (ke - epsilon));

				// Atoms need to escape, so what needs to happen is, we need to find which cell in 
				// All Atom Agents array is full, with atom a1, or a0, and set those spaces equal to null. 

					if( a0Partner1 == a1)
					{	((Atom[]) a0.allatomAgents[idx])[0] = null;
						//((Atom[]) a1.allatomAgents[idx])[0] = null;
					}  // Did I change the class field ?

					if( a0Partner2 == a1)
					{	((Atom[]) a0.allatomAgents[idx])[1] = null; }  // Did I change the class field ?

					if( a1Partner1 == a0)
					{	((Atom[]) a1.allatomAgents[idx])[0] = null; }  // Did I change the class field ?

					if( a1Partner2 == a1)
					{	((Atom[]) a1.allatomAgents[idx])[1] = null; }  // Did I change the class field ?


					nudge = eps;

					//System.out.println("escape");

/*					
					if(a0Partner1 == a1)
					{	
						(a0.allatomAgents[idx]) = null;	
					
if (( (Atom[]) a0.allatomAgents[idx])[0] == a1)
 ((Atom[]) a0.allatomAgents[idx])[0] = null;
					{
					a0.allatomAgents[idx] = null;
					a1.allatomAgents[idx] = null;

*/
				
				
				}

			}//end if(2*r2...
        }
		else { 	//not bonded to each other
			//well collision; decide whether to bond or have hard repulsion
			if (a1Saturated || a0Saturated) { //at least one is already bound;
											  // repel
				lastCollisionVirial = reduced_m * bij;
				nudge = eps;
				//System.out.println("repel");
			} else { //neither is taken; bond to each other
				lastCollisionVirial = 0.5
						* reduced_m
						* (bij + Math.sqrt(bij * bij + 4.0 * r2 * epsilon
								/ reduced_m));
				if (a0Partner1 == null)
					((Atom[]) a0.allatomAgents[idx])[0] = a1;
				else
					((Atom[]) a0.allatomAgents[idx])[1] = a1;
				if (a1Partner1 == null)
					((Atom[]) a1.allatomAgents[idx])[0] = a0;
				else
					((Atom[]) a1.allatomAgents[idx])[1] = a0;

				nudge = -eps;
				//System.out.println("bonded");
			}

		} //end if(a1Partner == a2) else

		lastCollisionVirialr2 = lastCollisionVirial / r2;
		dv.Ea1Tv1(lastCollisionVirialr2, dr);
		((ICoordinateKinetic) pair.atom0.coord).velocity().PEa1Tv1( rm0, dv);
		((ICoordinateKinetic) pair.atom1.coord).velocity().PEa1Tv1(-rm1, dv);
		a0.coord.position().PEa1Tv1(-falseTime * rm0, dv);
		a1.coord.position().PEa1Tv1(falseTime * rm1, dv);
		
		if (nudge != 0) 
		{
			a0.coord.position().PEa1Tv1(-nudge, dr);
			a1.coord.position().PEa1Tv1(nudge, dr);
		}

	}//end of bump
}//end of class
