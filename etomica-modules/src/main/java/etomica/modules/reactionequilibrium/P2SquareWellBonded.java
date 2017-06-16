/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.reactionequilibrium;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.atom.AtomLeafAgentManager;
import etomica.potential.P2SquareWell;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Energy;


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

    private static final long serialVersionUID = 1L;
	private double barrier;
    protected AtomLeafAgentManager agentManager;

	// *** Below are different types of constructor functions that Make P2SqaredWells ***	

//	public P2SquareWellBonded(Space space, double coreDiameter, double lambda, double epsilon) {
//		super(space, coreDiameter, lambda, epsilon);
		//constructs 2 BondChangeData classes, one for each atom in possible
		// bondchange event
//	}

	public P2SquareWellBonded(Space space, AtomLeafAgentManager aam, double coreDiameter, double lambda, double epsilon, boolean ignoreOverlap) {
		super(space, coreDiameter, lambda, epsilon, ignoreOverlap);
        agentManager = aam;
	}

	public void setBarrier(double b) {
		barrier = b;
	}

	public double getBarrier() {
		return barrier;
	}

	public Dimension getBarrierDimension() {
		return Energy.DIMENSION;
	}

	/**
	 * Computes next time of collision of the two atoms, assuming free-flight
	 * kinematics.
	 */

    public double collisionTime(IAtomList atoms, double falseTime) {
        
        if (ignoreOverlap) {
            
            IAtomKinetic atom0 = (IAtomKinetic)atoms.getAtom(0);
            IAtomKinetic atom1 = (IAtomKinetic)atoms.getAtom(1);

            // ** Makes 2 things, and atomPair pair, 
            IAtom a0Partner = (IAtom)agentManager.getAgent(atom0);
            if (a0Partner != atom1) {

                dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());

                dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
                dr.PEa1Tv1(falseTime,dv);
                boundary.nearestImage(dr);

                double r2 = dr.squared();

                if (r2 < wellDiameterSquared) {
                    //inside well but not mutually bonded; collide now if approaching, never if not.
                    double bij = dr.dot(dv);
                    return bij < 0.0 ? falseTime : Double.POSITIVE_INFINITY;
                }
            }
        }
        //mutually bonded, or outside well; collide as SW
        double time = super.collisionTime(atoms, falseTime);
//      if(!Double.isInfinite(time)) System.out.println("Collision time: "+time+" for "+atoms.toString());
        return time;
    }

	
	public void bump(IAtomList pair, double falseTime) {
		
		// *** Data Declaration Section

        IAtom atom0 = pair.getAtom(0);
        IAtom atom1 = pair.getAtom(1);
        IAtomKinetic coord0 = (IAtomKinetic)atom0;
        IAtomKinetic coord1 = (IAtomKinetic)atom1;
        dv.Ev1Mv2(coord1.getVelocity(), coord0.getVelocity());
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

		double r2 = dr.squared();
		double bij = dr.dot(dv);
		double nudge = 0;
		double eps = 1.0e-10;
		
		// ke is kinetic energy due to components of velocity
        double rm0 = atom0.getType().rm();
        double rm1 = atom1.getType().rm();
//		System.out.println("Bumping "+pair.toString());
		double reduced_m = 2.0 / (rm0 + rm1);
		double ke = bij * bij * reduced_m / (4.0 * r2);
		
		IAtom a0Partner = (IAtom)agentManager.getAgent(atom0);
		IAtom a1Partner = (IAtom)agentManager.getAgent(atom1);

		boolean a0Saturated = (a0Partner != null);
		boolean a1Saturated = (a1Partner != null);
		lastEnergyChange = 0.0;

		if (a0Partner == atom1) {	//atoms are bonded to each
            if (a1Partner != atom0) {
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
				lastEnergyChange = epsilon;

				// Atoms need to escape, so what needs to happen is, we need to find which cell in 
				// All Atom Agents array is full, with atom a1, or a0, and set those spaces equal to null. 

					if(a0Partner == atom1)
					{	
					    agentManager.setAgent(atom0, null);
                        agentManager.setAgent(atom1, null);
					}

					nudge = eps;
				
				}

			}//end if(2*r2...
        }
		else { 	//not bonded to each other
			//well collision; decide whether to bond or have hard repulsion
			if (a1Saturated || a0Saturated || lambda <= 1) { //at least one is already bound;
											  // repel
				lastCollisionVirial = reduced_m * bij;
				nudge = eps;
				//System.out.println("repel");
			} else { //neither is taken; bond to each other
				lastCollisionVirial = 0.5
						* reduced_m
						* (bij + Math.sqrt(bij * bij + 4.0 * r2 * epsilon
								/ reduced_m));
                lastEnergyChange = -epsilon;
                agentManager.setAgent(atom0, atom1);
                agentManager.setAgent(atom1, atom0);
				nudge = -eps;
				//System.out.println("bonded");
			}

		} //end if(a1Partner == a2) else

		lastCollisionVirialr2 = lastCollisionVirial / r2;
		dv.Ea1Tv1(lastCollisionVirialr2, dr);
		coord0.getVelocity().PEa1Tv1( rm0, dv);
		coord1.getVelocity().PEa1Tv1(-rm1, dv);
		coord0.getPosition().PEa1Tv1(-falseTime * rm0, dv);
		coord1.getPosition().PEa1Tv1(falseTime * rm1, dv);
		
		if (nudge != 0) 
		{
			coord0.getPosition().PEa1Tv1(-nudge, dr);
			coord1.getPosition().PEa1Tv1(nudge, dr);
		}

	}//end of bump
    
}//end of class
