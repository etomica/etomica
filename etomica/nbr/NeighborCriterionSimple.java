package etomica.nbr;

import etomica.Atom;
import etomica.Debug;
import etomica.NearestImageTransformerVector;
import etomica.NearestImageVectorSource;
import etomica.Phase;
import etomica.Space;
import etomica.nbr.cell.AtomsetIteratorCellular;
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
	
	public void setNeighborRadius(double r) {
		if (r < interactionRange) throw new IllegalArgumentException("Neighbor radius must be larger than interaction range");
		neighborRadius2 = r*r;
		double displacementLimit = (r - interactionRange) * safetyFactor;
		displacementLimit2 = displacementLimit * displacementLimit;
        r2MaxSafe = displacementLimit2 / (4.0*safetyFactor*safetyFactor);
	}
	
	public double getNeighborRadius() {
		return Math.sqrt(neighborRadius2);
	}
	
	public boolean needUpdate(Atom atom) {
		r2 = atom.coord.position().Mv1Squared((Vector)atom.allatomAgents[agentIndex]);
        if (Debug.ON && Debug.DEBUG_NOW && Debug.thisAtom(atom)) {
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
    
    public void setCellIterator(AtomsetIteratorCellular api) {
        api.getNbrCellIterator().setPeriod(phase.boundary().dimensions());
    }

    public void setNearestImageVectorSource(NearestImageVectorSource nivs) {
        nearestImageVectorSource = nivs;
    }
    
	public boolean unsafe() {
		if (Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
			System.out.println("some atom exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor));
		}
		return r2 > r2MaxSafe;
	}

	public boolean accept(Atom[] a) {
        nearestImageTransformer.setNearestImageVector(nearestImageVectorSource.getNearestImageVector());
		cPair.reset(a[0].coord,a[1].coord);
		if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.anyAtom(a)) || (Debug.LEVEL == 1 && Debug.allAtoms(a)))) {
			if (cPair.r2() < neighborRadius2 || (Debug.LEVEL > 1 && Debug.allAtoms(a))) {
				System.out.println("Atom "+a[0]+" and "+a[1]+" are "+(cPair.r2() < neighborRadius2 ? "" : "not ")+"neighbors, r2="+cPair.r2());
            }
		}
		return cPair.r2() < neighborRadius2;
	}
	
	public void reset(Atom atom) {
		((Vector)atom.allatomAgents[agentIndex]).E(atom.coord.position());
	}

	private double interactionRange, displacementLimit2, neighborRadius2;
	CoordinatePair cPair;
	protected static int agentIndex = Atom.requestAgentIndex(new Atom.AgentSource() {
		public Object makeAgent(Atom atom) {
			return atom.coord.position().clone();
		}
	});
	protected double safetyFactor;
	protected double r2, r2MaxSafe;
    private Phase phase;
    private final NearestImageTransformerVector nearestImageTransformer;
    private NearestImageVectorSource nearestImageVectorSource;
	
}
