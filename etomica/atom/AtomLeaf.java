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
