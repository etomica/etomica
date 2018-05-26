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
    private final int[] latticeIndex;

    public NeighborIteratorSite(NeighborSiteManager siteManager, Box box) {
        this.siteManager = siteManager;
        this.lattice = siteManager.getLattice();
        this.realIterator = new RectangularLatticeNbrIteratorAdjacent(lattice.D());
        realIterator.setLattice(lattice);
        realIterator.setPeriodicity(box.getBoundary().getPeriodicity());
        this.latticeIndex = new int[box.getSpace().getD()];
    }

    @Override
    public void forEachNeighbor(IAtom targetAtom, IteratorDirective.Direction direction, AtomPairConsumer upAction, AtomPairConsumer downAction) {
        lattice.latticeIndex(siteManager.getSite(targetAtom).latticeArrayIndex, latticeIndex);
        realIterator.setSite(this.latticeIndex);

        if (direction != IteratorDirective.Direction.DOWN) {
            realIterator.setDirection(IteratorDirective.Direction.UP);
            realIterator.reset();

            while (this.realIterator.hasNext()) {
                AtomSite site = ((AtomSite) realIterator.next());
                upAction.accept(targetAtom, site.getAtom());
            }
        }

        if (direction != IteratorDirective.Direction.UP) {
            realIterator.setDirection(IteratorDirective.Direction.DOWN);
            realIterator.reset();

            while (this.realIterator.hasNext()) {
                AtomSite site = ((AtomSite) realIterator.next());
                downAction.accept(site.getAtom(), targetAtom);
            }

        }
    }
}
