package etomica.nbr.cell;

import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomSequencerFactory;
import etomica.nbr.site.PotentialMasterSite;
import etomica.space.Space;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */

/*
 * Created on May 27, 2005
 */
public class PotentialMasterCell extends PotentialMasterSite {

    /**
     * Constructs with null AtomPositionDefinition, which indicates the position
     * definition given with each atom's AtomType should be used.
     * 
     * @param space the governing Space
     * @param range the neighbor distance.  May be changed after construction.
     */
    public PotentialMasterCell(Space space, double range) {
        this(space, range, null);
    }

    /**
     * @param space
     * @param positionDefinition
     */
    public PotentialMasterCell(Space space, double range,
            AtomPositionDefinition positionDefinition) {
        super(space, positionDefinition, new Api1ACell(space.D(), range));
    }

    public void setRange(double d) {
        ((Api1ACell)neighborIterator).getNbrCellIterator().setRange(d);
    }
    
    public AtomSequencerFactory sequencerFactory() {return AtomSequencerCell.FACTORY;}
    

}
