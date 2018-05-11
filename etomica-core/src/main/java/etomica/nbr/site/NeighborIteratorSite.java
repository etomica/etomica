package etomica.nbr.site;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLatticeNbrIteratorAdjacent;
import etomica.nbr.NeighborIterator;
import etomica.potential.IteratorDirective;

public class NeighborIteratorSite implements NeighborIterator {
    private final NeighborSiteManager siteManager;
    private final CellLattice lattice;
    private final RectangularLatticeNbrIteratorAdjacent realIterator;

    public NeighborIteratorSite(NeighborSiteManager siteManager, Box box) {
        this.siteManager = siteManager;
        this.lattice = siteManager.getLattice();
        this.realIterator = new RectangularLatticeNbrIteratorAdjacent(lattice.D());
        realIterator.setLattice(lattice);
        realIterator.setPeriodicity(box.getBoundary().getPeriodicity());
    }

    @Override
    public void forEachNeighbor(IAtom targetAtom, IteratorDirective.Direction direction, AtomPairConsumer upAction, AtomPairConsumer downAction) {
        // TODO
    }
}
