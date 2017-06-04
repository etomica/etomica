/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.molecule;

import etomica.api.*;
import etomica.atom.IMoleculePositionDefinition;
import etomica.atom.MoleculePositionGeometricCenter;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.atom.MoleculeAgentManager;
import etomica.atom.MoleculeAgentManager.MoleculeAgentSource;
import etomica.box.BoxAgentManager;
import etomica.box.BoxAgentSourceMoleculeManager;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Length;

/**
 * Simple neighbor criterion based on distance between the two molecules 
 * If the molecule moves further than r^2; the update method will kick in
 * to update the neighborlist.
 * 
 * 
 * @author Tai Boon Tan
 *
 */
public class CriterionSimpleMolecular implements NeighborCriterionMolecular, MoleculeAgentSource {

	public CriterionSimpleMolecular(Simulation sim, Space _space, double interactionRange, double neighborRadius) {
		super();
        this.space = _space;
        dr = space.makeVector();
		this.interactionRange = interactionRange;
        neighborRadius2 = neighborRadius * neighborRadius;
        setSafetyFactor(0.4);
        boxAgentManager = new BoxAgentManager<MoleculeAgentManager>(new BoxAgentSourceMoleculeManager(this, sim),MoleculeAgentManager.class, sim);
        moleculeSite = new MoleculePositionGeometricCenter(_space);
	}
	
	public void setMoleculePosition(IMoleculePositionDefinition positionDefinition){
		moleculeSite = positionDefinition; 
	}
	
	public void setSafetyFactor(double f) {
		if (f <= 0.0 || f >= 0.5) throw new IllegalArgumentException("safety factor must be positive and less than 0.5");
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
	
	public boolean needUpdate(IMolecule molecule) {
//        if (Debug.ON && interactionRange > Math.sqrt(neighborRadius2)) {
//            throw new IllegalStateException("Interaction range ("+interactionRange+") must be less than neighborRange ("+Math.sqrt(neighborRadius2)+")");
//        }
        
        
		r2 = moleculeSite.position(molecule).Mv1Squared((Vector)agentManager.getAgent(molecule));
		
//        if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.allAtoms(new AtomSetSinglet(atom))) {
//            System.out.println("atom "+atom+" displacement "+r2+" "+atom.getPosition());
//        }
//        
//        
//		if (Debug.ON && Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
//			System.out.println("atom "+atom+" exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor)+")");
//			System.out.println("old position "+agentManager.getAgent(atom));
//			System.out.println("new position "+atom.getPosition());
////            throw new RuntimeException("stop that");
//		}
		return r2 > displacementLimit2;
	
	
	}

	public void setBox(Box box) {
        boundary = box.getBoundary();
        agentManager = boxAgentManager.getAgent(box);
	}
    
	public boolean unsafe() {
//		if (Debug.ON && Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
//			System.out.println("some atom exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor));
//		}
		return r2 > r2MaxSafe;
	}

	public boolean accept(IMoleculeList pair) {
		dr.E(moleculeSite.position(pair.getMolecule(1)));
		dr.ME(moleculeSite.position(pair.getMolecule(0)));

		boundary.nearestImage(dr);
//        if (Debug.ON && neighborRadius2 < interactionRange*interactionRange) {
//            throw new IllegalStateException("neighbor radius "+Math.sqrt(neighborRadius2)+" is less than interaction range "+interactionRange);
//        }
//		if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.anyAtom(pair)) || (Debug.LEVEL == 1 && Debug.allAtoms(pair)))) {
//            double r2l = dr.squared(); 
//			if (r2l < neighborRadius2 || (Debug.LEVEL > 1 && Debug.allAtoms(pair))) {
//				System.out.println("Atom "+pair.getAtom(0)+" and "+pair.getAtom(1)+" are "+(r2l < neighborRadius2 ? "" : "not ")+"neighbors, r2="+r2l);
//            }
//		}
       return dr.squared() < neighborRadius2;
	}
	
	public void reset(IMolecule molecule) {
		((Vector)agentManager.getAgent(molecule)).E(moleculeSite.position(molecule));
	}

    public Class getMoleculeAgentClass() {
        return dr.getClass();
    }
    
    public Object makeAgent(IMolecule molecule) {
        return space.makeVector();
    }
    
    public void releaseAgent(Object agent, IMolecule molecule) {}

    protected final Space space;
    protected double interactionRange, displacementLimit2, neighborRadius2;
	protected final Vector dr;
    protected IBoundary boundary;
	protected double safetyFactor;
	protected double r2, r2MaxSafe;
    protected MoleculeAgentManager agentManager;
    protected final BoxAgentManager<MoleculeAgentManager> boxAgentManager;
    protected IMoleculePositionDefinition moleculeSite;

}
