/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.potential.P2HardGeneric;
import etomica.space.Vector;
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
public class P2SquareWellRadicalFasterer extends P2HardGeneric {

    protected final AtomLeafAgentManager<IAtom[]> agentManager;
    protected final IRandom random;
    protected double combinationProbability;
    protected double solventThermoFrac;

    public P2SquareWellRadicalFasterer(AtomLeafAgentManager<IAtom[]> aam,
                                       double coreDiameter, double lambda, double epsilon, IRandom random) {
        super(new double[]{coreDiameter, coreDiameter * lambda}, new double[]{Double.POSITIVE_INFINITY, -epsilon}, true);
        agentManager = aam;
        setSolventThermoFrac(0);
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
        IAtom[] nbrs = agentManager.getAgent(a);
        for (int i = 0; i < nbrs.length - 1; ++i) {
            if (nbrs[i] == null) {
                return false;
            }
        }
        return nbrs[nbrs.length - 1] == null;
    }

    protected boolean isEmpty(IAtom a) {
        return agentManager.getAgent(a)[0] == null;
    }

    /**
     * This will tell you what the lowest open space is in atom a
     */
    protected int lowest(IAtom a) {
        IAtom[] nbrs = agentManager.getAgent(a);
        int j = nbrs.length;
        for (int i = 0; i != j; ++i) {
            if (nbrs[i] == null) {
                return i;
            }
        }
        return j;
    }

    /**
     * This function tells you if two atoms are bonded
     */
    protected boolean areBonded(IAtom a, IAtom b) {
        IAtom[] nbrs = agentManager.getAgent(a);
        int j = nbrs.length;
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
        int i = lowest(a);
        int j = lowest(b);
        agentManager.getAgent(a)[i] = b;
        agentManager.getAgent(b)[j] = a;
    }

    /**
     * this function will makes a and b unreactive by setting both to be bonded
     * to themselves (a-a, b-b).
     */
    protected void disproportionate(IAtom a, IAtom b) {
        int i = lowest(a);
        int j = lowest(b);
        agentManager.getAgent(a)[i] = a;
        agentManager.getAgent(b)[j] = b;
    }


    public double collisionTime(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12, Vector v12, int collisionState) {

        double bij = r12.dot(v12);

        if (collisionState < 0) collisionState = 2;
        double r2 = r12.squared();
        if (bij < 0.0 && r2 < collisionDistances2[1] && !areBonded(atom1, atom2)) {
            // overlapped collide now
            return 0.001 * Math.sqrt(r2 / v12.squared());
        }
        return super.collisionTime(atom1, atom2, r12, v12, collisionState);
    }

    @Override
    protected int decideBump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, boolean core, double ke, double reducedMass, double bij, double r2, double[] du, double[] virial, double falseTime) {
        int newState = oldState + (core ? -1 : +1);
        if (areBonded(atom1, atom2)) {
            //atoms are bonded to each other -- there is no escape (Mu Ha Ha Ha!)
            virial[0] = 2.0 * reducedMass * bij;
            newState = oldState;
        } else {
            //not bonded to each other
            //well collision; decide whether to a) bond b) hard core repulsion
            // c) disproportionate
            boolean radical1 = isRadical(atom1);
            boolean radical2 = isRadical(atom2);
            boolean empty1 = isEmpty(atom1);
            boolean empty2 = isEmpty(atom2);
            double uJump = getEnergyForState(newState) - getEnergyForState(oldState);
            if (!radical1 && !radical2) {
                virial[0] = 2.0 * reducedMass * bij;
                newState = oldState;
            } else if (radical1 && radical2) {
                //radcial + radical.  terminate
                // if we're here and empty, we're initiator radical.  at least one of the atoms
                // has to be monomer radical
                if ((!empty1 || !empty2) && random.nextDouble() < combinationProbability) {
                    // react (bond)
                    virial[0] = reducedMass * (bij + (core ? +1 : -1) * Math.sqrt(bij * bij - 2.0 * r2 * uJump * solventThermoFrac / reducedMass));
                    bond(atom1, atom2);
                } else {
                    // disproportionate -- atoms are no longer reactive
                    virial[0] = 2.0 * reducedMass * bij;
                    newState = oldState;
                    disproportionate(atom1, atom2);
                }
            } else if ((radical1 && empty2) || (radical2 && empty1)) {
                //one is a radical, the other is a monomer.
                virial[0] = reducedMass * (bij + (core ? +1 : -1) * Math.sqrt(bij * bij - 2.0 * r2 * uJump * solventThermoFrac / reducedMass));
                bond(atom1, atom2);
            } else {
                // one of them is full
                newState = oldState;
                virial[0] = 2.0 * reducedMass * bij;
            }
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
}
