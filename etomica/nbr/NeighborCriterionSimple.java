package etomica.nbr;

import etomica.Atom;
import etomica.Debug;
import etomica.Phase;
import etomica.Space;

/**
 * Simple neighbor criterion based on distance moved by a leaf atom since
 * the last update.
 * @author andrew
 *
 */
public class NeighborCriterionSimple extends NeighborCriterion  {

	public NeighborCriterionSimple(Space space, double interactionRange, double neighborRadius) {
		super();
		this.interactionRange = interactionRange;
		double displacementLimit = (neighborRadius - interactionRange) * safetyFactor;
		displacementLimit2 = displacementLimit * displacementLimit;
		neighborRadius2 = neighborRadius * neighborRadius;
        cPair = space.makeCoordinatePair();
	}
	
	public void setSafetyFactor(double f) {
		if (safetyFactor < 0.0 || safetyFactor>0.5) throw new IllegalArgumentException("safety factor must be positive and less than 0.5");
		safetyFactor = f;
		double displacementLimit = (Math.sqrt(neighborRadius2) - interactionRange) * f;
		displacementLimit2 = displacementLimit * displacementLimit;
	}
	
	public double getSafetyFactor() {
		return safetyFactor;
	}
	
	public void setNeighborRadius(double r) {
		if (r < interactionRange) throw new IllegalArgumentException("Neighbor radius must be larger than interaction range");
		neighborRadius2 = r*r;
		double displacementLimit = (r - interactionRange) * safetyFactor;
		displacementLimit2 = displacementLimit * displacementLimit;
	}
	
	public double getNeighborRadius() {
		return Math.sqrt(neighborRadius2);
	}
	
	public boolean needUpdate(Atom atom) {
		r2 = atom.coord.position().Mv1Squared((Space.Vector)atom.allatomAgents[agentIndex]);
		if (Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
			System.out.println("atom "+atom+" exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor));
			System.out.println("old position "+(Space.Vector)atom.allatomAgents[agentIndex]);
			System.out.println("new position "+atom.coord.position());
		}
		return r2 > displacementLimit2;
	}

	public void setPhase(Phase phase) {
    	cPair.setBoundary(phase.boundary());
	}

	public boolean unsafe() {
		if (Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
			System.out.println("some atom exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor));
		}
		return r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor);
	}

	public boolean accept(Atom[] a) {
		cPair.reset(a[0].coord,a[1].coord);
		if (Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.anyAtom(a)) || (Debug.LEVEL == 1 && Debug.allAtoms(a)))) {
			if (cPair.r2() < neighborRadius2 || (Debug.LEVEL > 1 && Debug.allAtoms(a))) 
				System.out.println("Atom "+a[0]+" and "+a[1]+" are "+(cPair.r2() < neighborRadius2 ? "" : "not ")+"neighbors, r2="+cPair.r2());
		}
		return cPair.r2() < neighborRadius2;
	}
	
	public void reset(Atom atom) {
		((Space.Vector)atom.allatomAgents[agentIndex]).E(atom.coord.position());
	}

	private double interactionRange, displacementLimit2, neighborRadius2;
	private Space.CoordinatePair cPair;
	protected static int agentIndex = Atom.requestAgentIndex(new Atom.AgentSource() {
		public Object makeAgent(Atom atom) {
			return atom.coord.position().clone();
		}
	});
	protected double safetyFactor = 0.4;
	protected double r2;
	
}
