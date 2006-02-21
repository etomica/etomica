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
import etomica.space.CoordinatePair;
import etomica.space.Vector;
import etomica.util.Debug;

/**
 * Simple neighbor criterion based on distance moved by a leaf atom since
 * the last update.
 * @author andrew
 *
 */
public class CriterionSimple extends NeighborCriterion implements AgentSource  {

	public CriterionSimple(Simulation sim, double interactionRange, double neighborRadius) {
		super();
		this.interactionRange = interactionRange;
        neighborRadius2 = neighborRadius * neighborRadius;
        setSafetyFactor(0.4);
        cPair = new CoordinatePair(sim.space);
        phaseAgentManager = new PhaseAgentManager(new PhaseAgentSourceAtomManager(this),sim.speciesRoot);
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
	
	public boolean isRangeDependent() {
		return true;
	}
	
	public boolean needUpdate(Atom atom) {
		r2 = ((AtomLeaf)atom).coord.position().Mv1Squared(agents[atom.getGlobalIndex()]);
        if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.allAtoms(atom)) {
            System.out.println("atom "+atom+" displacement "+r2+" "+((AtomLeaf)atom).coord.position());
        }
		if (Debug.ON && Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
			System.out.println("atom "+atom+" exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor)+")");
			System.out.println("old position "+agents[atom.getGlobalIndex()]);
			System.out.println("new position "+((AtomLeaf)atom).coord.position());
            throw new RuntimeException("stop that");
		}
		return r2 > displacementLimit2;
	}

	public void setPhase(Phase phase) {
        cPair.setNearestImageTransformer(phase.getBoundary());
        agentManager = (AtomAgentManager[])phaseAgentManager.getAgents();
        agents = (Vector[])agentManager[phase.getIndex()].getAgents();
	}
    
	public boolean unsafe() {
		if (Debug.ON && Debug.DEBUG_NOW && r2 > displacementLimit2 / (4.0*safetyFactor*safetyFactor)) {
			System.out.println("some atom exceeded safe limit ("+r2+" > "+displacementLimit2 / (4.0*safetyFactor*safetyFactor));
		}
		return r2 > r2MaxSafe;
	}

	public boolean accept(AtomSet pair) {
		cPair.reset((AtomPair)pair);
		if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.anyAtom(pair)) || (Debug.LEVEL == 1 && Debug.allAtoms(pair)))) {
			if (cPair.r2() < neighborRadius2 || (Debug.LEVEL > 1 && Debug.allAtoms(pair))) {
				System.out.println("Atom "+((AtomPair)pair).atom0+" and "+((AtomPair)pair).atom1+" are "+(cPair.r2() < neighborRadius2 ? "" : "not ")+"neighbors, r2="+cPair.r2());
            }
		}
		return cPair.r2() < neighborRadius2;
	}
	
	public void reset(Atom atom) {
		agents[atom.getGlobalIndex()].E(((AtomLeaf)atom).coord.position());
	}

    public Class getAgentClass() {
        return cPair.dr().getClass();
    }
    
    public Object makeAgent(Atom atom) {
        return atom.type.isLeaf() ? ((AtomLeaf)atom).coord.position().clone() : null;
    }
    
    public void releaseAgent(Object agent, Atom atom) {}

    private double interactionRange, displacementLimit2, neighborRadius2;
	private final CoordinatePair cPair;
	protected double safetyFactor;
	protected double r2, r2MaxSafe;
    private AtomAgentManager[] agentManager;
    protected Vector[] agents;
    private final PhaseAgentManager phaseAgentManager;
}
