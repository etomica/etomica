package etomica.nbr;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomSet;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.phase.PhaseAgentSourceAtomManager;
import etomica.simulation.Simulation;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.util.Debug;

/**
 * Simple neighbor criterion based on distance moved by a leaf atom since
 * the last update.  The potential is assumed to be from a wall that exists
 * at the phase boundaries or at some fixed position within the phase.
 * @author andrew
 */
public class CriterionPositionWall implements NeighborCriterion, AgentSource, java.io.Serializable {

	public CriterionPositionWall(Simulation sim) {
		super();
		this.interactionRange = Double.NaN;
        this.neighborRange = Double.NaN;
        setBoundaryWall(true);
        setSafetyFactor(0.8);
        phaseAgentManager = new PhaseAgentManager(new PhaseAgentSourceAtomManager(this),sim.speciesRoot);
	}

    /**
     * Sets the safety factor (between 0 and 1.0) that determines when the 
     * criterion thinks it needs an update.  Safety factors near 1.0 allow
     * atoms to travel farther before needing an update, but are more risky.
     */
	public void setSafetyFactor(double f) {
		if (f <= 0.0 || f >= 1.0) throw new IllegalArgumentException("safety factor must be positive and less than 1.0");
		safetyFactor = f;
        rMaxSafe = (neighborRange - interactionRange);
		displacementLimit = rMaxSafe * safetyFactor;
	}

    /**
     * returns the safety factor
     */
    public double getSafetyFactor() {
        return safetyFactor;
    }

    /**
     * Sets the orientation of the wall to be perpendicular to the given dimension
     */
    public void setWallDim(int d) {
        neighborDim = d;
    }
    
    /**
     * Returns the interaction range of the wall potential.
     */
    public double getInteractionRange() {
        return interactionRange;
    }
	
    /**
     * Informs the criterion of the interaction range of the wall potential.
     */
    public void setInteractionRange(double r) {
        interactionRange = r;
        if (neighborRange == Double.NaN) {
            return;
        }
        rMaxSafe = (neighborRange - interactionRange);
        displacementLimit = rMaxSafe * safetyFactor;
    }        

    public Dimension getInteractionRangeDimension() {
        return Length.DIMENSION;
    }
    
    /**
     * Sets the neighbor range of the criterion.  Atoms within the given 
     * distance of the wall are "accepted".
     */
	public void setNeighborRange(double r) {
		if (interactionRange != Double.NaN && r < interactionRange) throw new IllegalArgumentException("Neighbor radius must be larger than interaction range");
		neighborRange = r;
        if (interactionRange == Double.NaN) {
            return;
        }
        rMaxSafe = (neighborRange - interactionRange);
        displacementLimit = rMaxSafe * safetyFactor;
	}
    
    public double getNeighborRange() {
        return neighborRange;
    }
    
    public Dimension getNeighborRangeDimension() {
        return Length.DIMENSION;
    }
    
	/**
     * Returns true if the walls are at the phase boundaries.
     */
    public boolean isBoundaryWall() {
        return isBoundaryWall;
    }

    /**
     * Sets whether the walls are at the phase boundaries or not.
     */
    public void setBoundaryWall(boolean isBoundaryWall) {
        this.isBoundaryWall = isBoundaryWall;
    }

    /**
     * Sets the position of the wall.  This parameter is ignored if 
     * isBoundaryWall is true.
     */
    public void setWallPosition(double p) {
        wallPosition = p;
    }
    
    /**
     * Returns the position of the wall.  This parameter is ignored if
     * isBoundaryWall is true.
     */
    public double getWallPosition() {
        return wallPosition;
    }
    
    public Dimension getWallPositionDimension() {
        return Length.DIMENSION;
    }

	public boolean needUpdate(Atom atom) {
        dr = Math.abs(((AtomLeaf)atom).coord.position().x(neighborDim) - ((DoubleWrapper)agentManager.getAgent(atom)).x);
        if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.allAtoms(atom)) {
            System.out.println("atom "+atom+" displacement "+dr+" "+((AtomLeaf)atom).coord.position());
        }
		if (Debug.ON && Debug.DEBUG_NOW && dr > rMaxSafe) {
			System.out.println("atom "+atom+" exceeded safe limit ("+dr+" > "+rMaxSafe+")");
			System.out.println("old position "+((DoubleWrapper)agentManager.getAgent(atom)).x);
			System.out.println("new position "+((AtomLeaf)atom).coord.position().x(neighborDim));
            throw new RuntimeException("stop that");
		}
		return dr > displacementLimit;
	}

	public void setPhase(Phase phase) {
        boxSize = phase.getBoundary().getDimensions().x(neighborDim);
        agentManager = (AtomAgentManager)phaseAgentManager.getAgent(phase);
	}
    
	public boolean unsafe() {
		if (Debug.ON && Debug.DEBUG_NOW && dr > rMaxSafe) {
			System.out.println("some atom exceeded safe limit ("+dr+" > "+rMaxSafe+")");
		}
		return dr > rMaxSafe;
	}

	public boolean accept(AtomSet atom) {
		dr = ((AtomLeaf)atom).coord.position().x(neighborDim);
        if (!isBoundaryWall) {
            dr = Math.abs(dr - wallPosition);
        }
        else {
            if (dr > 0.0) {
                dr = 0.5*boxSize - dr;
            }
            else {
                dr = dr + 0.5*boxSize;
            }
        }
		if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 0 && Debug.allAtoms(atom))) {
			if (dr < neighborRange || Debug.LEVEL > 1) {
				System.out.println("Atom "+atom+" is "+(dr < neighborRange ? "" : "not ")+"interacting, dr="+dr);
            }
		}
		return dr < neighborRange;
	}
	
	public void reset(Atom atom) {
		((DoubleWrapper)agentManager.getAgent(atom)).x = ((AtomLeaf)atom).coord.position().x(neighborDim);
	}

    public Class getAgentClass() {
        return DoubleWrapper.class;
    }
    
    public Object makeAgent(Atom atom) {
        return atom.type.isLeaf() ? new DoubleWrapper() : null;
    }
    
    public void releaseAgent(Object agent, Atom atom) {}

    protected static class DoubleWrapper implements java.io.Serializable {
        public double x;
    }
    
    private double interactionRange, displacementLimit, neighborRange;
    private int neighborDim;
    private boolean isBoundaryWall;
    private double wallPosition;
    private double boxSize;
	protected double safetyFactor;
	protected double dr, rMaxSafe;
    protected AtomAgentManager agentManager;
    private final PhaseAgentManager phaseAgentManager;
}
