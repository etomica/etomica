/*
 * Created on May 27, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.nbr.cell;

import etomica.Space;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomSequencerFactory;
import etomica.nbr.site.PotentialMasterSite;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class PotentialMasterCell extends PotentialMasterSite {

    /**
     * @param space
     */
    public PotentialMasterCell(Space space) {
        this(space,null);
    }

    /**
     * @param space
     * @param positionDefinition
     */
    public PotentialMasterCell(Space space,
            AtomPositionDefinition positionDefinition) {
        super(space, positionDefinition, new Api1ACell(space.D()));
    }

    public void setRange(double d) {
        ((Api1ACell)neighborIterator).getNbrCellIterator().setRange(d);
    }
    
    public AtomSequencerFactory sequencerFactory() {return AtomSequencerCell.FACTORY;}
    

}
