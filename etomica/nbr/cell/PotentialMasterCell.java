package etomica.nbr.cell;

import etomica.atom.Atom;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.site.PotentialMasterSite;
import etomica.phase.Phase;
import etomica.potential.Potential;
import etomica.potential.PotentialCalculation;
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
     * Creates PotentialMasterCell with default (1.0) range.  Range
     * should be set manually via setRange method.
     */
    public PotentialMasterCell(Space space) {
        this(space,1.0);
    }
    
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
        super.setRange(d);
        ((Api1ACell)neighborIterator).getNbrCellIterator().setNeighborDistance(d);
    }
    
    public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
        currentNeighborCellManager = getNbrCellManager(phase);
        super.calculate(phase,id,pc);
    }
    protected void calculate(Atom atom, IteratorDirective id, PotentialCalculation pc, final Potential[] potentials) {
        ((Api1ACell)neighborIterator).setCentralCell(currentNeighborCellManager.getCell(atom));
        super.calculate(atom,id,pc,potentials);
    }
    public AtomSequencerFactory sequencerFactory() {return AtomSequencerCell.FACTORY;}
    
    private NeighborCellManager currentNeighborCellManager;
}
