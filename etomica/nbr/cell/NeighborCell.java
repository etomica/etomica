/*
 * History
 * Created on Nov 23, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.Atom;
import etomica.AtomFactory;
import etomica.AtomList;
import etomica.AtomSequencerSimple;
import etomica.AtomTreeNodeGroup;
import etomica.AtomTreeNodeLeaf;
import etomica.AtomType;
import etomica.Space;
import etomica.lattice.Primitive;
import etomica.lattice.Site;

/**
 * @author kofke
 *
 * Site used to form array of cells for cell-based neighbor listing.  Each
 * cell is capable of holding lists of atoms that are in them.
 */
public class NeighborCell extends Site {

    //invoked by factory class
    private NeighborCell(Space space, AtomType type, AtomTreeNodeGroup parent) {
        super(space, type, parent);
    }
    
    private AtomList[] occupants;

    /**
     * Returns a factory that can make NeighborCell instances.
     * @param space space in which cells reside
     * @param primitive provides polytope that defines size and shape of each cell
     * @return a new NeighborCell factory
     */
    public static AtomFactory makeFactory(Space space, Primitive primitive) {
        return new NeighborCellFactory(space, primitive);
    }
    private static class NeighborCellFactory extends AtomFactory {
        
        public NeighborCellFactory(Space space, Primitive primitive) {
            super(space, AtomSequencerSimple.FACTORY, AtomTreeNodeLeaf.FACTORY);
            atomType = new AtomTypeCell(this, primitive.unitCell());
        }
        public boolean isGroupFactory() {return false;}
        protected Atom build(AtomTreeNodeGroup parent) {
            return new NeighborCell(space, atomType, parent);
        }
        public Atom build(Atom atom) {
            if(!(atom instanceof NeighborCell)) throw new IllegalArgumentException("Illegal atom for build");
            return atom;
        }
        
    }
}
