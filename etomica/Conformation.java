package etomica;

import etomica.atom.AtomList;

/**
 * General class for assignment of coordinates to a group of atoms.
 */
 
 /* History of changes
  * 09/01/02 (DAK) added field to flag whether total momentum should be set to
  * zero when initalizing (made in conjunction with change to Space.
  * CoordinateGroup classes, which no longer do randomizeMomentum to zero total
  * momentum). 
  * 09/04/02 (DAK) added Null inner class and NULL field
  * 01/21/04 (DAK) added initializeCoordinate(Phase) method
  */
public abstract class Conformation {

    public Conformation(Space space) {
        this.space = space;
    }

    public abstract void initializePositions(AtomList atomList);

    protected final Space space;
}//end of Configuration
