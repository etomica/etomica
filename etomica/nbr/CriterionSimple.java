package etomica.nbr;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.phase.PhaseAgentSourceAtomManager;
import etomica.simulation.Simulation;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.util.Debug;

/**
 * Simple neighbor criterion based on distance moved by a leaf atom since
 * the last update.
 * @author andrew
 *
 */
public class CriterionSimple implements NeighborCriterion, AgentSource, java.io.Serializable {

	public CriterionSimple(Simulation sim, double interactionRange, double neighborRadius) {
		super();
        space = sim.getSpace();
        dr = space.makeVector();
		this.interactionRange = interactionRange;
        neighborRadius2 = neighborRadius * neighborRadius;
        setSafetyFactor(0.4);
        phaseAgentManager = new PhaseAgentManager(new PhaseAgentSourceAtomManager(this),sim.getSpeciesRoot());
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
	
	public boolean needUpdate(Atom atom) {
        if (Debug.ON && interactionRange > Math.sqrt(neighborRadius2)) {
            throw new IllegalStateException("Interaction range ("+interactionRange+") must be less than neighborRange ("+Math.sqrt(neighborRadius2)+")");
        }
		r2 = ((AtomLeaf)atom).getCoord().position().Mv1Squared((Vector)agentManager.getAgent(atom));
        if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.allAtoms(atom)) {
            System.out.println("atom "+atom+" displacement "+r2+" "+((AtomLeaf)atom).getCoord().position());
        }
		if (Debug.ON && Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
			System.out.println("atom "+atom+" exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor)+")");
			System.out.println("old position "+agentManager.getAgent(atom));
			System.out.println("new position "+((AtomLeaf)atom).getCoord().position());
            throw new RuntimeException("stop that");
		}
		return r2 > displacementLimit2;
	}

	public void setPhase(Phase phase) {
        nearestImageTransformer = phase.getBoundary();
        agentManager = (AtomAgentManager)phaseAgentManager.getAgent(phase);
	}
    
	public boolean unsafe() {
		if (Debug.ON && Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
			System.out.println("some atom exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor));
		}
		return r2 > r2MaxSafe;
	}

	public boolean accept(AtomSet pair) {
        dr.Ev1Mv2(((AtomLeaf)((AtomPair)pair).atom1).getCoord().position(),((AtomLeaf)((AtomPair)pair).atom0).getCoord().position());
        nearestImageTransformer.nearestImage(dr);
        if (Debug.ON && neighborRadius2 < interactionRange*interactionRange) {
            throw new IllegalStateException("neighbor radius "+Math.sqrt(neighborRadius2)+" is less than interaction range "+interactionRange);
        }
		if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.anyAtom(pair)) || (Debug.LEVEL == 1 && Debug.allAtoms(pair)))) {
            double r2l = dr.squared(); 
			if (r2l < neighborRadius2 || (Debug.LEVEL > 1 && Debug.allAtoms(pair))) {
				System.out.println("Atom "+((AtomPair)pair).atom0+" and "+((AtomPair)pair).atom1+" are "+(r2l < neighborRadius2 ? "" : "not ")+"neighbors, r2="+r2l);
            }
		}
		return dr.squared() < neighborRadius2;
	}
	
	public void reset(Atom atom) {
        ((Vector)agentManager.getAgent(atom)).E(((AtomLeaf)atom).getCoord().position());
	}

    public Class getAgentClass() {
        return dr.getClass();
    }
    
    public Object makeAgent(Atom atom) {
        return atom.getType().isLeaf() ? space.makeVector() : null;
    }
    
    public void releaseAgent(Object agent, Atom atom) {}

    private static final long serialVersionUID = 1L;
    protected final Space space;
    private double interactionRange, displacementLimit2, neighborRadius2;
	private final Vector dr;
    private NearestImageTransformer nearestImageTransformer;
	protected double safetyFactor;
	protected double r2, r2MaxSafe;
    private AtomAgentManager agentManager;
    private final PhaseAgentManager phaseAgentManager;
}
