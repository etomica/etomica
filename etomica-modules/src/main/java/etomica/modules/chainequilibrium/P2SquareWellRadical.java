/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.potential.P2SquareWell;
import etomica.space.Space;
import etomica.util.random.IRandom;


/**
 * Similar to square-well potential, but considers and alters bonding states
 * with collisions based on free radical reactions.  Fully reacted atoms are
 * inert.  A radical can react (bond) with a monomer.  If two radicals meet,
 * they can bond or disproportionate (become unreactive without forming a bond)
 * based on the combinationProbability.
 * 
 * @author Andrew Schultz
 */
public class P2SquareWellRadical extends P2SquareWell {

    private static final long serialVersionUID = 1L;
    protected final AtomLeafAgentManager agentManager;
    protected final IRandom random;
    protected double combinationProbability;

    public P2SquareWellRadical(Space space, AtomLeafAgentManager aam,
                               double coreDiameter, double lambda, double epsilon, IRandom random) {
        super(space, coreDiameter, lambda, epsilon, true);
        agentManager = aam;
        this.random = random;
    }

    /**
     * Sets the probability that two radicals will bond instead of becoming
     * unreactive.
     */
    public void setCombinationProbability(double newCombinationProbability) {
        if (newCombinationProbability < 0 || newCombinationProbability > 1) {
            throw new IllegalArgumentException("invalid probability");
        }
        combinationProbability = newCombinationProbability;
    }

    /**
     * Returns the probability that two radicals will bond instead of becoming
     * unreactive.
     */
    public double getCombinationProbability() {
        return combinationProbability;
    }

    /**
     * This function return true if the given atom is a radical (meaning that
     * it has only one unreacted site).
     */
    protected boolean isRadical(IAtom a) {
        IAtom[] nbrs = (IAtom[])agentManager.getAgent(a);
        for(int i=0; i < nbrs.length-1; ++i){
            if (nbrs[i] == null) {
                return false;
            }
        }
        return nbrs[nbrs.length-1] == null;
    }

    protected boolean isEmpty(IAtom a) {
        return ((IAtom[])agentManager.getAgent(a))[0] == null;
    }

    /**
     * This will tell you what the lowest open space is in atom a
     */
    protected int lowest(IAtom a){
        IAtom[] nbrs = (IAtom[])agentManager.getAgent(a);
        int j = nbrs.length;
        for(int i=0; i != j; ++i){
            if (nbrs[i] == null) {
                return i;
            }
        }
        return j; 
    }

    /**
     * This function tells you if two atoms are bonded
     */
    protected boolean areBonded(IAtom a, IAtom b){
        IAtom[] nbrs = (IAtom[])agentManager.getAgent(a);
        int j = nbrs.length;
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
        int i = lowest(a);
        int j = lowest(b);
        ((IAtom[])agentManager.getAgent(a))[i] = b;
        ((IAtom[])agentManager.getAgent(b))[j] = a;
    }

    /**
     * this function will makes a and b unreactive by setting both to be bonded
     * to themselves (a-a, b-b).
     */
    protected void disproportionate(IAtom a, IAtom b){
        int i = lowest(a);
        int j = lowest(b);
        ((IAtom[])agentManager.getAgent(a))[i] = a;
        ((IAtom[])agentManager.getAgent(b))[j] = b;
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
            if (r2 < wellDiameterSquared) {
                boolean areBonded = areBonded(atom0, atom1);
                if (!areBonded) {
                    //inside well but not mutually bonded; collide now if approaching
                    return (dr.dot(dv) < 0) ? falseTime : Double.POSITIVE_INFINITY;
                }
            }
        }
        //mutually bonded, or outside well; collide as SW
        return super.collisionTime(atoms, falseTime);
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

        double rm0 = atom0.getType().rm();
        double rm1 = atom1.getType().rm();
        
        double reduced_m = 2.0 /  + (rm0 + rm1);

        IAtom atomLeaf0 = atom0;
        IAtom atomLeaf1 = atom1;
        if (areBonded(atomLeaf0,atomLeaf1)) {		//atoms are bonded to each other
            lastCollisionVirial = reduced_m * bij;
            if (2 * r2 > (coreDiameterSquared + wellDiameterSquared)) {															
                // there is no escape (Mu Ha Ha Ha!), nudge back inside
                nudge = -eps;
            }
        }
        else { 	//not bonded to each other
            //well collision; decide whether to a) bond b) hard core repulsion
            // c) disproportionate
            boolean radical0 = isRadical(atomLeaf0);
            boolean radical1 = isRadical(atomLeaf1);
            boolean empty0 = isEmpty(atomLeaf0);
            boolean empty1 = isEmpty(atomLeaf1);
            if (!radical0 && !radical1) {
                lastCollisionVirial = reduced_m * bij;
                nudge = eps;
            }
            else if (radical0 && radical1) {
                //radcial + radical.  terminate
                // if we're here and empty, we're initiator radical.  at least one of the atoms
                // has to be monomer radical
                if ((!empty0 || !empty1) && random.nextDouble() < combinationProbability) {
                    // react (bond)
                    lastCollisionVirial = 0.5* reduced_m* (bij + Math.sqrt(bij * bij + 4.0 * r2 *2* epsilon/ reduced_m));
                    bond(atomLeaf0, atomLeaf1);
                    nudge = -eps;
			    }
                else {
                    // disproportionate
                    lastCollisionVirial = reduced_m * bij;
                    disproportionate(atomLeaf0, atomLeaf1);
                    nudge = eps;
                }
            }
            else if ((radical0 && empty1) || (radical1 && empty0)) {
                //one is a radical, the other is a monomer.
                lastCollisionVirial = 0.5* reduced_m* (bij + Math.sqrt(bij * bij + 4.0 * r2 * epsilon/ reduced_m));
                bond(atomLeaf0, atomLeaf1);
                nudge = -eps;
            }
            else {
                // one of them is full
                lastCollisionVirial = reduced_m * bij;
                nudge = eps;
            }
        }

        lastCollisionVirialr2 = lastCollisionVirial / r2;
        dv.Ea1Tv1(lastCollisionVirialr2, dr);
        atom0.getVelocity().PEa1Tv1(rm0, dv);
        atom1.getVelocity().PEa1Tv1(-rm1, dv);
        atom0.getPosition().PEa1Tv1(-falseTime * rm0, dv);
        atom1.getPosition().PEa1Tv1(falseTime * rm1, dv);

        if (nudge != 0) {
            atom0.getPosition().PEa1Tv1(-nudge, dr);
            atom1.getPosition().PEa1Tv1(nudge, dr);
        }
    }
}
