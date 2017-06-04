/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.api.IBoundary;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomSetSinglet;
import etomica.box.BoxAgentManager;
import etomica.box.BoxAgentSourceAtomManager;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.util.Debug;

/**
 * Simple neighbor criterion based on distance moved by a leaf atom since
 * the last update.
 * @author andrew
 *
 */
public class CriterionSimple implements NeighborCriterion, AgentSource<Vector> {

	public CriterionSimple(Simulation sim, Space _space, double interactionRange, double neighborRadius) {
		super();
        this.space = _space;
        dr = space.makeVector();
		this.interactionRange = interactionRange;
        neighborRadius2 = neighborRadius * neighborRadius;
        setSafetyFactor(0.4);
        BoxAgentSourceAtomManager<Vector> basam = new BoxAgentSourceAtomManager<Vector>(this, Vector.class);
        boxAgentManager = new BoxAgentManager<AtomLeafAgentManager<Vector>>(basam,AtomLeafAgentManager.class,sim);
	}
	
	public void setSafetyFactor(double f) {
		if (f <= 0.0 || f > 0.5) throw new IllegalArgumentException("safety factor must be positive and less than 0.5");
		safetyFactor = f;
		double displacementLimit = (Math.sqrt(neighborRadius2) - interactionRange) * f;
		displacementLimit2 = displacementLimit * displacementLimit;
        r2MaxSafe = displacementLimit2 / (4.0*safetyFactor*safetyFactor);
	}
	
	public double getSafetyFactor() {
		return safetyFactor;
	}
	
	public void setNeighborRange(double r) {
		neighborRadius2 = r*r;
		double displacementLimit = (r - interactionRange) * safetyFactor;
		displacementLimit2 = displacementLimit * displacementLimit;
        r2MaxSafe = displacementLimit2 / (4.0*safetyFactor*safetyFactor);
	}
    
    public double getNeighborRange() {
        return Math.sqrt(neighborRadius2);
    }
    
    public Dimension getNeighborRangeDimension() {
        return Length.DIMENSION;
    }
    
    public void setInteractionRange(double newInteractionRange) {
        interactionRange = newInteractionRange;
        double displacementLimit = (Math.sqrt(neighborRadius2) - interactionRange) * safetyFactor;
        displacementLimit2 = displacementLimit * displacementLimit;
        r2MaxSafe = displacementLimit2 / (4.0*safetyFactor*safetyFactor);
    }
    
    public double getInteractionRange() {
        return interactionRange;
    }
    
    public Dimension getInteractionRangeDimension() {
        return Length.DIMENSION;
    }
	
	public boolean needUpdate(IAtom atom) {
        if (Debug.ON && interactionRange > Math.sqrt(neighborRadius2)) {
            throw new IllegalStateException("Interaction range ("+interactionRange+") must be less than neighborRange ("+Math.sqrt(neighborRadius2)+")");
        }
		r2 = atom.getPosition().Mv1Squared(agentManager.getAgent(atom));
        if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.allAtoms(new AtomSetSinglet(atom))) {
            System.out.println("atom "+atom+" displacement "+r2+" "+atom.getPosition());
        }
		if (Debug.ON && Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
			System.out.println("atom "+atom+" exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor)+")");
			System.out.println("old position "+agentManager.getAgent(atom));
			System.out.println("new position "+atom.getPosition());
//            throw new RuntimeException("stop that");
		}
		return r2 > displacementLimit2;
	}

	public void setBox(Box box) {
        boundary = box.getBoundary();
        agentManager = boxAgentManager.getAgent(box);
	}
    
	public boolean unsafe() {
		if (Debug.ON && Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
			System.out.println("some atom exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor));
		}
		return r2 > r2MaxSafe;
	}

	public boolean accept(IAtomList pair) {
        dr.Ev1Mv2(pair.getAtom(1).getPosition(),pair.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        if (Debug.ON && neighborRadius2 < interactionRange*interactionRange) {
            throw new IllegalStateException("neighbor radius "+Math.sqrt(neighborRadius2)+" is less than interaction range "+interactionRange);
        }
		if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.anyAtom(pair)) || (Debug.LEVEL == 1 && Debug.allAtoms(pair)))) {
            double r2l = dr.squared(); 
			if (r2l < neighborRadius2 || (Debug.LEVEL > 1 && Debug.allAtoms(pair))) {
				System.out.println("Atom "+pair.getAtom(0)+" and "+pair.getAtom(1)+" are "+(r2l < neighborRadius2 ? "" : "not ")+"neighbors, r2="+r2l);
            }
		}
		return dr.squared() < neighborRadius2;
	}
	
	public void reset(IAtom atom) {
        agentManager.getAgent(atom).E(atom.getPosition());
	}

    public Vector makeAgent(IAtom atom, Box agentBox) {
        Vector v = space.makeVector();
        // atom isn't necessarily in the position.  but if atom-adding code is smart,
        // it will be in the appropriate position.
        v.E(atom.getPosition());
        return v;
    }
    
    public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {}

    protected final Space space;
    protected double interactionRange, displacementLimit2, neighborRadius2;
	protected final Vector dr;
    protected IBoundary boundary;
	protected double safetyFactor;
	protected double r2, r2MaxSafe;
    protected AtomLeafAgentManager<Vector> agentManager;
    protected final BoxAgentManager<AtomLeafAgentManager<Vector>> boxAgentManager;
}
