/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.potential.P2HardGeneric;
import etomica.potential.Potential2Soft;
import etomica.space.Vector;


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
public class P2SquareWellBonded extends P2HardGeneric implements Potential2Soft {

    protected final AtomLeafAgentManager<IAtom[]> agentManager;
    protected double solventThermoFrac;

    public P2SquareWellBonded(AtomLeafAgentManager<IAtom[]> aam, double coreDiameter, double lambda, double epsilon) {
        super(new double[]{coreDiameter, coreDiameter * lambda}, new double[]{Double.POSITIVE_INFINITY, -epsilon}, true);
        agentManager = aam;
        setSolventThermoFrac(0);
        ringResult = new RingResult();
    }

    /**
     * This function will tell the user, if passed an atom whether or not that atom can bond
     */
    protected boolean full(IAtom a) {
        IAtom[] nbrs = agentManager.getAgent(a);
        int j = nbrs.length;    //check INDEXING
        for (int i = 0; i != j; ++i) {
            if (nbrs[i] == null) {
                return false;
            }
        }
        return true;
    }

    /**
     * This will tell you what the lowest open space is in atom a
     */
    protected int lowest(IAtom a) {
        IAtom[] nbrs = agentManager.getAgent(a);
        int j = nbrs.length;    //check INDEXING
        for (int i = 0; i != j; ++i) {
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
    protected boolean areBonded(IAtom a, IAtom b) {
        IAtom[] nbrs = agentManager.getAgent(a);
        int j = nbrs.length;    //check INDEXING
        for (int i = 0; i != j; ++i) {
            if (nbrs[i] == b) {
                return true;
            }
        }
        return false;
    }

    /**
     * this function will bond atoms a & b together
     */
    protected void bond(IAtom a, IAtom b) {
        if (areBonded(a, b)) {            // Error Checking, what about double bonds?
            throw new RuntimeException(a + " and " + b + " are already bonded");
        }
        int i = lowest(a);        // (0 is the First Space)
        int j = lowest(b);
        agentManager.getAgent(a)[i] = b;
        agentManager.getAgent(b)[j] = a;
    }

    /**
     * this function will bond atoms a & b together
     */
    protected void checkRing(IAtom a, IAtom b, int maxBondCount) {
        IAtom[] aNbrs = agentManager.getAgent(a);
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
            IAtom[] nextNbrs = agentManager.getAgent(next);
            if (nextNbrs.length == 3) {
                // encountered a cross-linker.  rings are OK.
                ringResult.linker = next;
                ringResult.bondCount = bondCount;
                return;
            } else if (nextNbrs.length == 1) {
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
    protected void unbond(IAtom a, IAtom b) {
        if (!areBonded(a, b)) {        // Error Checking
            throw new RuntimeException(a + " and " + b + " are not bonded");
        }
        boolean success = false;
        // Unbonding the Atom, Atom A's side
        IAtom[] nbrs = agentManager.getAgent(a);
        for (int i = 0; i < nbrs.length; ++i) {
            if (nbrs[i] == b) {
                nbrs[i] = null;
                success = true;
            }
        }
        if (!success) {
            throw new RuntimeException("oops #1 " + b + " not in " + a + " list");
        }
        success = false;
        // Unbonding the Atom, Atom B's side
        nbrs = agentManager.getAgent(b);
        for (int i = 0; i < nbrs.length; ++i) {
            if (nbrs[i] == a) {
                nbrs[i] = null;
                success = true;
            }
        }
        if (!success) {
            throw new RuntimeException("oops #2 " + b + " not in " + a + " list");
        }
    }

    @Override
    protected int decideBump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, boolean core, double ke, double reducedMass, double bij, double r2, double[] du, double[] virial, double falseTime) {
        int newState = oldState + (core ? -1 : +1);
        boolean bonded = areBonded(atom1, atom2);
        boolean canBond = !bonded && bij < 0 && !full(atom1) && !full(atom2);
        if (canBond) {
            int maxRingBonds = 20;
            IAtom[] aNbrs = agentManager.getAgent(atom1);
            if (aNbrs.length == 3) {
                // cross linker
                checkRing(atom2, atom1, maxRingBonds);
                canBond = ringResult.linker == atom1 || (ringResult.linker != null && ringResult.foundRing);
            } else {
                checkRing(atom1, atom2, maxRingBonds);
//    		        System.out.println("checkRing "+atom0+" "+atom1+" linker0 "+ringResult.linker);
                //                System.out.println(atom0+" "+atom1+" "+ringBonds);
                if (ringResult.linker != null) {
                    IAtom linker0 = ringResult.linker;
                    int ringBonds0 = ringResult.bondCount;
                    checkRing(atom2, atom1, maxRingBonds - ringBonds0);
//    	                System.out.println("checkRing "+atom0+" "+atom1+" linker1 "+ringResult.linker);
                    if (ringResult.linker == linker0) {
                        // ring contains only one linker and is too small
                        canBond = false;
                    }
                } else {
                    canBond = !ringResult.foundRing;
                }
            }
        }
        double uJump = (canBond || bonded) ? (getEnergyForState(atom1, atom2, newState) - getEnergyForState(atom1, atom2, oldState)) : Double.POSITIVE_INFINITY;
        if (ke < uJump) {
            // not enough ke; bounce off core
            virial[0] = 2.0 * reducedMass * bij;
            newState = oldState;
        } else {
            if (bonded) {
                // capture
                unbond(atom1, atom2);
            } else {
                bond(atom1, atom2);
            }
            virial[0] = reducedMass * (bij + (core ? +1 : -1) * Math.sqrt(bij * bij - 2.0 * r2 * uJump * solventThermoFrac / reducedMass));
            du[0] = uJump;
        }
        return newState;
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

    @Override
    public Vector[][] gradientAndTorque(IAtomList atoms) {
        return new Vector[0][];
    }

    protected static class RingResult {
        public IAtom linker;
        public int bondCount;
        public boolean foundRing;
    }

    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        int s = getState(atom1, atom2, dr12);
        if (areBonded(atom1, atom2)) {
            return getEnergyForState(atom1, atom2, s);
        }
        return (0 <= s && s < 2) ? Double.POSITIVE_INFINITY : 0;
    }
}

