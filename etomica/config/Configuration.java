package etomica.config;

import etomica.Phase;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.space.Space;

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

    public Configuration(Space space) {
        this.space = space;
    }

    public abstract void initializePositions(AtomList[] atomList);
    
    public void initializeCoordinates(Phase phase) {
    	setDimensions(phase.boundary().dimensions().toArray());
		AtomList speciesAgentList = ((AtomTreeNodeGroup)phase.getSpeciesMaster().node).childList;
        AtomIteratorListSimple speciesAgentIterator = new AtomIteratorListSimple(speciesAgentList);
        AtomList[] moleculeLists = new AtomList[speciesAgentList.size()];
        int i=0;
        speciesAgentIterator.reset();
        while (speciesAgentIterator.hasNext()) {
            SpeciesAgent agent = (SpeciesAgent)speciesAgentIterator.nextAtom();
            moleculeLists[i++] = ((AtomTreeNodeGroup)agent.node).childList;
        }
        initializePositions(moleculeLists);
    }
    
    public void setDimensions(double[] dim) {
        dimensions = (double[])dim.clone();
    }

    public double[] getDimensions() {return dimensions;}
    
    protected double[] dimensions;
    protected Space space;
    
}//end of Configuration
