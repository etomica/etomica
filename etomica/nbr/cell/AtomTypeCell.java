/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.AtomFactory;
import etomica.AtomType;
import etomica.math.geometry.Polytope;

/**
 * @author kofke
 *
 * AtomType used by the cells of the cell neighbor-list lattice.
 */
public class AtomTypeCell extends AtomType {

    /**
     * @param creator
     */
    public AtomTypeCell(AtomFactory creator, Polytope cell) {
        super(creator);
        this.unitCell = cell;
     }

    public final Polytope unitCell;
}
