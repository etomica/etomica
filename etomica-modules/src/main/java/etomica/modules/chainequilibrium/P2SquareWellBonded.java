/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.AtomLeafAgentManager;
import etomica.potential.P2SquareWell;
import etomica.space.Space;


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
    protected final AtomLeafAgentManager agentManager;
    protected Box box;
    protected double solventThermoFrac;

	public P2SquareWellBonded(Space space, AtomLeafAgentManager aam, double coreDiameter, double lambda, double epsilon) {
		super(space, coreDiameter, lambda, epsilon, true);
        agentManager = aam;
        setSolventThermoFrac(0);
        ringResult = new RingResult();
	}

    public void setBox(Box newBox){
        box = newBox;
        super.setBox(box);
    }
    
	/**
     * This function will tell the user, if passed an atom weither or not that atom can bond
	 */
	protected boolean full(IAtom a) {
        IAtom[] nbrs = (IAtom[])agentManager.getAgent(a);
		int j = nbrs.length;	//check INDEXING
		for(int i=0; i != j; ++i){
			if (nbrs[i] == null) {
				return false;
			}
		}
		return true; 
	}
	
	/**
     * This will tell you what the lowest open space is in atom a
	 */
	protected int lowest(IAtom a){
        IAtom[] nbrs = (IAtom[])agentManager.getAgent(a);
		int j = nbrs.length;	//check INDEXING
		for(int i=0; i != j; ++i){
			if (nbrs[i] == null) {
				return i;
			}
		}
		return j; 
	}
	
	/**
     * This function tells you if two atoms are bonded
     * This could probably be public, although a public version would
     * need to first re-retrieve agents
	 */
	protected boolean areBonded(IAtom a, IAtom b){
        IAtom[] nbrs = (IAtom[])agentManager.getAgent(a);
		int j = nbrs.length;	//check INDEXING
		for(int i=0; i != j; ++i){
			if (nbrs[i] == b){		
				return true;
			}
		}
		return false; 	
	}
    
	/**
     * this function will bond atoms a & b together
	 */
	protected void bond(IAtom a, IAtom b){
		if (areBonded(a,b)){			// Error Checking, what about double bonds?
			throw new RuntimeException(a+" and "+b+" are already bonded");
		}
		int i = lowest(a);		// (0 is the First Space) 
		int j = lowest(b);
        ((IAtom[])agentManager.getAgent(a))[i] = b;
        ((IAtom[])agentManager.getAgent(b))[j] = a;
	}
	
    /**
     * this function will bond atoms a & b together
     */
    protected void checkRing(IAtom a, IAtom b, int maxBondCount){
        IAtom[] aNbrs = ((IAtom[])agentManager.getAgent(a));
        if (aNbrs.length < 2) {
            ringResult.linker = null;
            ringResult.foundRing = false;
            return;
        }
        IAtom next = aNbrs[0];
        if (next == null) {
            next = aNbrs[1];
        }
        if (next == null) {
            // a is unbonded
            ringResult.linker = null;
            ringResult.foundRing = false;
            return;
        }
        int bondCount = 1;
        IAtom prev = a;
        while (true) {
            IAtom[] nextNbrs = ((IAtom[])agentManager.getAgent(next));
            if (nextNbrs.length == 3) {
                // encountered a cross-linker.  rings are OK.
                ringResult.linker = next;
                ringResult.bondCount = bondCount;
                return;
            }
            else if (nextNbrs.length == 1) {
                // encountered a monofunctional (terminal) group.  so, no ring.
                ringResult.linker = null;
                ringResult.foundRing = false;
                return;
            }
            IAtom nextNext = nextNbrs[0];
            if (nextNext == prev) {
                // we want |next|'s bonded partner's bonded partner that isn't |next|
                nextNext = nextNbrs[1];
            }
            if (nextNext == null) {
                // termination.  no ring
                ringResult.linker = null;
                ringResult.foundRing = false;
                return;
            }
            bondCount++;
            if (bondCount > maxBondCount) {
                // might be a ring, but if it is, it's large enough to be OK
                ringResult.linker = null;
                ringResult.foundRing = false;
                return;
            }
            prev = next;
            next = nextNext;
            if (next == b) {
                // found a ring
                ringResult.linker = null;
                ringResult.foundRing = true;
                ringResult.bondCount = bondCount;
                return;
            }
        }
    }
    
	/**
     * this function unbonds two atoms
	 */
	protected void unbond(IAtom a, IAtom b){
		if (!areBonded(a,b)){		// Error Checking
            throw new RuntimeException(a+" and "+b+" are not bonded");
		}
        boolean success = false;
		// Unbonding the Atom, Atom A's side
        IAtom[] nbrs = (IAtom[])agentManager.getAgent(a);
		for(int i=0; i < nbrs.length; ++i){
			if (nbrs[i] == b){
				nbrs[i] = null;
                success = true;
			}
		}
        if (!success) {
            throw new RuntimeException("oops #1 "+b+" not in "+a+" list");
        }
        success = false;
		// Unbonding the Atom, Atom B's side
        nbrs = (IAtom[])agentManager.getAgent(b);
		for(int i=0; i < nbrs.length; ++i){
			if (nbrs[i] == a){
				nbrs[i] = null;
                success = true;
			}
		}
        if (!success) {
            throw new RuntimeException("oops #2 "+b+" not in "+a+" list");
        }
	}
	
	/**
	 * Computes next time of collision of the two atoms, assuming free-flight
	 * kinematics.
	 */
	public double collisionTime(IAtomList atoms, double falseTime) {
	
		if (ignoreOverlap) {
			
            IAtomKinetic atom0 = (IAtomKinetic)atoms.get(0);
            IAtomKinetic atom1 = (IAtomKinetic)atoms.get(1);
            dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
            
            dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
            dr.PEa1Tv1(falseTime,dv);
            boundary.nearestImage(dr);

			double r2 = dr.squared();
			double bij = dr.dot(dv);
			boolean areBonded = areBonded(atom0, atom1);
			//inside well but not mutually bonded; collide now if approaching
            if (!areBonded && r2 < wellDiameterSquared) {
                return (bij < 0) ? falseTime : Double.POSITIVE_INFINITY;
            }
		}
		//mutually bonded, or outside well; collide as SW
		double time = super.collisionTime(atoms, falseTime);
//		if(!Double.isInfinite(time)) System.out.println("Collision time: "+time+" for "+atoms.toString());
		return time;
	}

	
	public void bump(IAtomList pair, double falseTime) {

        IAtomKinetic atom0 = (IAtomKinetic)pair.get(0);
        IAtomKinetic atom1 = (IAtomKinetic)pair.get(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

		double r2 = dr.squared();
		double bij = dr.dot(dv);
		double nudge = 0;
		double eps = 1.0e-10;
		
		// ke is kinetic energy due to components of velocity
		
		double rm0 = atom0.getType().rm();
		double rm1 = atom1.getType().rm();
		
		double reduced_m = 2.0 /  + (rm0 + rm1);
		double ke = bij * bij * reduced_m / (4.0 * r2);
		
		IAtom atomLeaf0 = atom0;
        IAtom atomLeaf1 = atom1;
		if (areBonded(atomLeaf0,atomLeaf1)) {		//atoms are bonded to each
			if (2 * r2 < (coreDiameterSquared + wellDiameterSquared)) { // Hard-core collision															
				lastCollisionVirial = reduced_m * bij;
	
			} else { 				// Well collision assume separating because mutually bonded
				if (ke < epsilon) 		// Not enough kinetic energy to escape
				{ 
					lastCollisionVirial = reduced_m * bij;
					nudge = -eps;
				} 
				else{ 	
				    lastCollisionVirial = 0.5 * reduced_m * bij- Math.sqrt(reduced_m * r2 * (ke - epsilon*solventThermoFrac));
					unbond(atomLeaf1,atomLeaf0);
					nudge = eps;
				}

			}
        }
		else { 	//not bonded to each other
			//well collision; decide whether to bond or have hard repulsion
		    boolean canBond = !full(atomLeaf0) && !full(atomLeaf1);
		    if (canBond) {
                int maxRingBonds = 20;
		        IAtom[] aNbrs = ((IAtom[])agentManager.getAgent(atomLeaf0));
		        if (aNbrs.length == 3) {
		            // cross linker
                    checkRing(atomLeaf1, atomLeaf0, maxRingBonds);
                    canBond = ringResult.linker == atomLeaf0 || (ringResult.linker != null && ringResult.foundRing);
		        }
		        else {
    		        checkRing(atomLeaf0, atomLeaf1, maxRingBonds);
//    		        System.out.println("checkRing "+atom0+" "+atom1+" linker0 "+ringResult.linker);
    //                System.out.println(atom0+" "+atom1+" "+ringBonds);
    		        if (ringResult.linker != null) {
    		            IAtom linker0 = ringResult.linker;
    		            int ringBonds0 = ringResult.bondCount;
    		            checkRing(atomLeaf1, atomLeaf0, maxRingBonds - ringBonds0);
//    	                System.out.println("checkRing "+atom0+" "+atom1+" linker1 "+ringResult.linker);
    		            if (ringResult.linker == linker0) {
    		                // ring contains only one linker and is too small
    		                canBond = false;
    		            }
    		        }
    		        else {
    		            canBond = !ringResult.foundRing;
    		        }
		        }
		    }
			if (!canBond) { 
				lastCollisionVirial = reduced_m * bij;
				nudge = eps;
			}
			else {
			    //neither is taken; bond to each other
                lastCollisionVirial = 0.5* reduced_m* (bij + Math.sqrt(bij * bij + 4.0 * r2 * epsilon*solventThermoFrac/ reduced_m));
				bond(atom0, atom1);
				nudge = -eps;
			}
		} 

		lastCollisionVirialr2 = lastCollisionVirial / r2;
		dv.Ea1Tv1(lastCollisionVirialr2, dr);
		atom0.getVelocity().PEa1Tv1(rm0, dv);
		atom1.getVelocity().PEa1Tv1(-rm1, dv);
		atom0.getPosition().PEa1Tv1(-falseTime * rm0, dv);
		atom1.getPosition().PEa1Tv1(falseTime * rm1, dv);
		
		if (nudge != 0) 
		{
			atom0.getPosition().PEa1Tv1(-nudge, dr);
			atom1.getPosition().PEa1Tv1(nudge, dr);
		}
	}

    /**
     * Returns the fraction of well energy that is gained or lost by an atom
     * pair when they hop in or leave their well.
     */
    public double getSolventThermoFrac() {
        return 1 - solventThermoFrac;
    }

    /**
     * Sets the fraction of well energy that is gained by an atom pair when
     * they hop in their well.  This is also the fraction of energy they give
     * up when they leave the well.  The pair does still need to have
     * sufficient energy to escape the well at its full strength.
     */
    public void setSolventThermoFrac(double newSolventThermoFrac) {
        if (newSolventThermoFrac < 0 || newSolventThermoFrac > 1) {
            throw new IllegalArgumentException("0 <= value <= 1");
        }
        solventThermoFrac = 1 - newSolventThermoFrac;
    }

    protected final RingResult ringResult;
    
    protected static class RingResult {
        public IAtom linker;
        public int bondCount;
        public boolean foundRing;
    }
}

