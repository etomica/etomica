package etomica.nbr;

import etomica.Atom;
import etomica.AtomPair;
import etomica.Debug;
import etomica.NearestImageTransformerVector;
import etomica.NearestImageVectorSource;
import etomica.Phase;
import etomica.Space;
import etomica.atom.AtomPairVector;
import etomica.space.CoordinatePair;
import etomica.space.Vector;

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
        neighborRadius2 = neighborRadius * neighborRadius;
        setSafetyFactor(0.4);
        cPair = space.makeCoordinatePair();
        nearestImageTransformer = new NearestImageTransformerVector();
        nearestImageTransformer.setPlus(false);
        cPair.setNearestImageTransformer(nearestImageTransformer);
	}
	
	public void setSafetyFactor(double f) {
		if (safetyFactor < 0.0 || safetyFactor > 0.5) throw new IllegalArgumentException("safety factor must be positive and less than 0.5");
		safetyFactor = f;
		double displacementLimit = (Math.sqrt(neighborRadius2) - interactionRange) * f;
		displacementLimit2 = displacementLimit * displacementLimit;
        r2MaxSafe = displacementLimit2 / (4.0*safetyFactor*safetyFactor);
	}
	
	public double getSafetyFactor() {
		return safetyFactor;
	}
	
	public void setNeighborRange(double r) {
		if (r < interactionRange) throw new IllegalArgumentException("Neighbor radius must be larger than interaction range");
		neighborRadius2 = r*r;
		double displacementLimit = (r - interactionRange) * safetyFactor;
		displacementLimit2 = displacementLimit * displacementLimit;
        r2MaxSafe = displacementLimit2 / (4.0*safetyFactor*safetyFactor);
	}
	
	public double getNeighborRange() {
		return Math.sqrt(neighborRadius2);
	}
	
	public boolean needUpdate(Atom atom) {
		r2 = atom.coord.position().Mv1Squared((Vector)atom.allatomAgents[agentIndex]);
        if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.thisAtom(atom)) {
            System.out.println("atom "+atom+" displacement "+r2+" "+atom.coord.position());
        }
		if (Debug.ON && Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
			System.out.println("atom "+atom+" exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor)+")");
			System.out.println("old position "+atom.allatomAgents[agentIndex]);
			System.out.println("new position "+atom.coord.position());
            throw new RuntimeException("stop that");
		}
		return r2 > displacementLimit2;
	}

	public void setPhase(Phase phase) {
	    this.phase = phase;
	}
    
	public boolean unsafe() {
		if (Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
			System.out.println("some atom exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor));
		}
		return r2 > r2MaxSafe;
	}

	public boolean accept(AtomPair pair) {
        nearestImageTransformer.setNearestImageVector(((AtomPairVector)pair).nearestImageVector);
		cPair.reset(pair);
		if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.anyAtom(pair)) || (Debug.LEVEL == 1 && Debug.allAtoms(pair)))) {
			if (cPair.r2() < neighborRadius2 || (Debug.LEVEL > 1 && Debug.allAtoms(pair))) {
				System.out.println("Atom "+(pair).atom0+" and "+(pair).atom1+" are "+(cPair.r2() < neighborRadius2 ? "" : "not ")+"neighbors, r2="+cPair.r2());
            }
		}
		return cPair.r2() < neighborRadius2;
	}
	
	public void reset(Atom atom) {
		((Vector)atom.allatomAgents[agentIndex]).E(atom.coord.position());
	}

	private double interactionRange, displacementLimit2, neighborRadius2;
	private final CoordinatePair cPair;
	protected static final int agentIndex = Atom.requestAgentIndex(new Atom.AgentSource() {
		public Object makeAgent(Atom atom) {
			return (atom.coord != null) ? atom.coord.position().clone() : null;
		}
	});
	protected double safetyFactor;
	protected double r2, r2MaxSafe;
    private Phase phase;
    private final NearestImageTransformerVector nearestImageTransformer;
    private NearestImageVectorSource nearestImageVectorSource;
	
}
