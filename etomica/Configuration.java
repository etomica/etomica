package etomica;

import etomica.action.AtomActionAccelerateTo;
import etomica.action.AtomActionRandomizeVelocity;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.space.Vector;

/**
 * General class for assignment of coordinates to a group of atoms.
 * 
 * @author David Kofke
 */
 
 /* History of changes
  * 09/01/02 (DAK) added field to flag whether total momentum should be set to
  * zero when initalizing (made in conjunction with change to Space.
  * CoordinateGroup classes, which no longer do randomizeMomentum to zero total
  * momentum). 
  * 09/04/02 (DAK) added Null inner class and NULL field
  * 01/21/04 (DAK) added initializeCoordinate(Phase) method
  */
public abstract class Configuration implements java.io.Serializable {

    protected boolean zeroTotalMomentum = false;
    protected double[] dimensions;
    protected AtomGroupAction atomActionRandomizeVelocity;
    protected AtomActionAccelerateTo atomActionAccelerateTo;
    protected Space space;
    protected Vector work;
    
    public Configuration(Space space) {
        this.space = space;
        atomActionRandomizeVelocity = new AtomGroupAction(new AtomActionRandomizeVelocity());
        atomActionAccelerateTo = new AtomActionAccelerateTo(space);
        work = space.makeVector();
    }

    public abstract void initializePositions(AtomList atomList);

    /** 
     * Sets flag indicating if total momentum should be zeroed after initializing
     * sub-atom momenta.  Default is true.
     */
    public void setZeroTotalMomentum(boolean b) {zeroTotalMomentum = b;}
    /** 
     * Flag indicating if total momentum should be zeroed after initializing
     * sub-atom momenta.  Default is true.
     */
    public boolean isZeroTotalMomentum() {return zeroTotalMomentum;}
 /**   
  * All atom velocities are set such that all have the same total momentum (corresponding to
  * the default value of temperature), with the direction at random
  */
    public void initializeMomenta(Atom atom) {
        initializeMomenta(atom, Default.TEMPERATURE);
    }
    public void initializeMomenta(Atom atom, double temperature) {
        ((AtomActionRandomizeVelocity)atomActionRandomizeVelocity.getAction()).setTemperature(temperature);
        atomActionRandomizeVelocity.actionPerformed(atom);
        if(zeroTotalMomentum) {
            work.E(0.0);
            atomActionAccelerateTo.setTargetVelocity(work);
            atomActionAccelerateTo.actionPerformed(atom);
        }
    }//end of initializeMomenta
    
    /**
     * Initializes positions and momenta of the atoms in the given atom group.
     */
    public void initializeCoordinates(Atom group) {
        if(group.node.isLeaf()) throw new IllegalArgumentException("Error in Configuration.initializeCoordinates:  Attempt to initialize child coordinates of leaf atom");

        initializePositions(((AtomTreeNodeGroup)group.node).childList);
        if (space.isKinetic()) {
            initializeMomenta(group);
        }
    }
}//end of Configuration
