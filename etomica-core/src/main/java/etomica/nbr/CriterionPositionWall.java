/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomSetSinglet;
import etomica.box.BoxAgentManager;
import etomica.box.BoxAgentSourceAtomManager;
import etomica.nbr.CriterionPositionWall.DoubleWrapper;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.util.Debug;

/**
 * Simple neighbor criterion based on distance moved by a leaf atom since
 * the last update.  The potential is assumed to be from a wall that exists
 * at the box boundaries or at some fixed position within the box.
 * @author andrew
 */
public class CriterionPositionWall implements NeighborCriterion, AgentSource<DoubleWrapper> {

	public CriterionPositionWall(Simulation sim) {
		super();
		this.interactionRange = Double.NaN;
        this.neighborRange = Double.NaN;
        setBoundaryWall(true);
        setSafetyFactor(0.8);
        boxAgentManager = new BoxAgentManager<AtomLeafAgentManager<DoubleWrapper>>(new BoxAgentSourceAtomManager<DoubleWrapper>(this,DoubleWrapper.class),AtomLeafAgentManager.class,sim);
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
     * Returns true if the walls are at the box boundaries.
     */
    public boolean isBoundaryWall() {
        return isBoundaryWall;
    }

    /**
     * Sets whether the walls are at the box boundaries or not.
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

	public boolean needUpdate(IAtom atom) {
        dr = Math.abs(atom.getPosition().getX(neighborDim) - agentManager.getAgent(atom).x);
        if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.allAtoms(new AtomSetSinglet(atom))) {
            System.out.println("atom "+atom+" displacement "+dr+" "+atom.getPosition());
        }
		if (Debug.ON && Debug.DEBUG_NOW && dr > rMaxSafe) {
			System.out.println("atom "+atom+" exceeded safe limit ("+dr+" > "+rMaxSafe+")");
			System.out.println("old position "+ agentManager.getAgent(atom).x);
			System.out.println("new position "+(atom).getPosition().getX(neighborDim));
            System.err.println("stop that");
		}
		return dr > displacementLimit;
	}

	public void setBox(Box box) {
        boxSize = box.getBoundary().getBoxSize().getX(neighborDim);
        agentManager = boxAgentManager.getAgent(box);
	}
    
	public boolean unsafe() {
		if (Debug.ON && Debug.DEBUG_NOW && dr > rMaxSafe) {
			System.out.println("some atom exceeded safe limit ("+dr+" > "+rMaxSafe+")");
		}
		return dr > rMaxSafe;
	}

	public boolean accept(IAtomList atom) {
		dr = (atom.getAtom(0)).getPosition().getX(neighborDim);
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
	
	public void reset(IAtom atom) {
		agentManager.getAgent(atom).x = atom.getPosition().getX(neighborDim);
	}

    public DoubleWrapper makeAgent(IAtom atom, Box agentBox) {
        return new DoubleWrapper();
    }
    
    public void releaseAgent(DoubleWrapper agent, IAtom atom, Box agentBox) {}

    protected static class DoubleWrapper {
        public double x;
    }

    private double interactionRange, displacementLimit, neighborRange;
    private int neighborDim;
    private boolean isBoundaryWall;
    private double wallPosition;
    private double boxSize;
	protected double safetyFactor;
	protected double dr, rMaxSafe;
    protected AtomLeafAgentManager<DoubleWrapper> agentManager;
    private final BoxAgentManager<AtomLeafAgentManager<DoubleWrapper>> boxAgentManager;
}
