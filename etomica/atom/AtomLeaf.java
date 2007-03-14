package etomica.atom;

import etomica.space.CoordinateFactorySphere;
import etomica.space.ICoordinate;
import etomica.space.Space;

 /**
  * Object corresponding to one physical atom or group of atoms. Each atom holds
  * the following publicly accessible fields:
  * <ul>
  * <li>a Coordinate instance (fieldname: coord) that is constructed by the
  * governing space class; the coordinate stores information about the state of
  * the atom -- usually its position and momentum, but other definitions are possible
  * <li>an AtomType instance (fieldname: type) that holds information this atom
  * has in common with other atoms made by the same factory
  * <li>an instance of AtomTreeNode (fieldname: node) that is used to place it
  * in the species hierarchy
  * <li>an instance of AtomLinker (fieldname: seq, for sequencer) that places
  * it in the linked list of children held by its parenttom holds an object
  * <li>an object (fieldname: ia, for integrator agent) that may be used to
  * store information needed by the integrator
  * <li>an array of objects (field name: allAtomAgents) that can be used to
  * store in each atom any object needed by another class; such a class must
  * implement Atom.AgentSource and request its object be stored in every atom by
  * invoking Atom.requestAgentIndex before any atoms are constructed. The
  * integer returned by this method will indicate the location in the
  * allAtomAgents array where the agent-source's object will be held.
  * </ul>
  * <p>
  * @author David Kofke, Andrew Schultz, and C. Daniel Barnes
  * 
  */
public class AtomLeaf extends Atom {

    public AtomLeaf(ICoordinate coord, AtomType type) {
        super(type);
        this.coord = coord;
    }
    
    /**
     * Makes a simple atom for the given space.  Coordinate is non-kinetic sphere;
     * node is for a leaf atom; type is a sphere with unit mass and unit size, 
     * unique to the new atom; depth is 0.
     */
    public AtomLeaf(Space space) {
        super();
        coord = new CoordinateFactorySphere(space,false).makeCoordinate();
    }
    
    /**
     * @return this atom's coordinate
     */
    public ICoordinate getCoord() {
        return coord;
    }

    public boolean isLeaf() {return true;}
    
    /**
     * Returns this node's atom.
     */
    public AtomLeaf firstLeafAtom() {return this;}
    
    /**
     * Returns this node's atom.
     */
    public AtomLeaf lastLeafAtom() {return this;}
    
    /**
     * Returns 1.
     */
    public int leafAtomCount() {return 1;}
    
    /**
     * Returns 0.
     */
    public int childAtomCount() {return 0;}

    public final void setLeafIndex(int newLeafIndex) {
        leafIndex = newLeafIndex;
    }
    
    public final int getLeafIndex() {
        return leafIndex;
    }
    
    /**
     * leafIndex is an index to the AtomArrayList of all leaf atoms in the phase.
     * List is maintained by the speciesMaster node.
     */
    private int leafIndex;

    private static final long serialVersionUID = 2L;
    private final ICoordinate coord;
    
}
