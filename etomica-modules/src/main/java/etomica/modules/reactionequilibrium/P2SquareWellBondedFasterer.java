/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.reactionequilibrium;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.potential.P2HardGeneric;
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
public class P2SquareWellBondedFasterer extends P2HardGeneric {

	protected AtomLeafAgentManager<IAtom> agentManager;

	public P2SquareWellBondedFasterer(AtomLeafAgentManager<IAtom> aam, double coreDiameter, double lambda, double epsilon) {
		super(new double[]{coreDiameter, coreDiameter * lambda}, new double[]{Double.POSITIVE_INFINITY, -epsilon}, true);
		agentManager = aam;
	}

	public void determineBonding(IAtom atom1, IAtom atom2, Vector rij) {
		double r2 = rij.squared();
		IAtom a1Partner = agentManager.getAgent(atom1);
		IAtom a2Partner = agentManager.getAgent(atom2);
		boolean a1Saturated = (a1Partner != null);
		boolean a2Saturated = (a2Partner != null);
		if (a1Saturated || a2Saturated || r2 > collisionDistances2[1]) return;
		// form bond
		agentManager.setAgent(atom1, atom2);
		agentManager.setAgent(atom2, atom1);
	}

	@Override
	protected double[] getCollisionDistances2(IAtom atom1, IAtom atom2) {
		IAtom a1Partner = agentManager.getAgent(atom1);
		IAtom a2Partner = agentManager.getAgent(atom2);
		boolean a1Saturated = (a1Partner != null);
		boolean a2Saturated = (a2Partner != null);
		if ((a1Saturated || a2Saturated) && a1Partner != atom2) return new double[]{collisionDistances2[1]};
		return collisionDistances2;
	}

	@Override
	protected int decideBump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, boolean core, double ke, double reducedMass, double bij, double r2, double[] du, double[] virial, double falseTime) {
		IAtom a1Partner = agentManager.getAgent(atom1);
		IAtom a2Partner = agentManager.getAgent(atom2);

		boolean a1Saturated = (a1Partner != null);
		boolean a2Saturated = (a2Partner != null);
		boolean bonded = a1Partner == atom2;

		int newState = oldState + (core ? -1 : +1);
		double uJump = getEnergyForState(atom1, atom2, newState) - getEnergyForState(atom1, atom2, oldState);
		if (ke < uJump || ((a1Saturated || a2Saturated) && !bonded)) {
			// not enough ke; bounce off core
			virial[0] = 2.0 * reducedMass * bij;
			newState = oldState;
		} else {
			// capture or escape
			virial[0] = reducedMass * (bij + (core ? +1 : -1) * Math.sqrt(bij * bij - 2.0 * r2 * uJump / reducedMass));
			du[0] = uJump;

			if (bonded) {
				// unbond
				agentManager.setAgent(atom1, null);
				agentManager.setAgent(atom2, null);
			} else {
				// form bond
				agentManager.setAgent(atom1, atom2);
				agentManager.setAgent(atom2, atom1);
			}
		}
		return newState;
	}
}
