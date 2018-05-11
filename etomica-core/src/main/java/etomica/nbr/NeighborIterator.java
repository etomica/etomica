package etomica.nbr;

import etomica.atom.IAtom;
import etomica.potential.IteratorDirective;

import java.util.function.BiConsumer;
import java.util.function.Consumer;

public interface NeighborIterator {

    void forEachNeighbor(IAtom targetAtom, IteratorDirective.Direction direction, AtomPairConsumer upAction, AtomPairConsumer downAction);

    default void forEachNeighbor(IAtom targetAtom, IteratorDirective.Direction direction, AtomPairConsumer action) {
        forEachNeighbor(targetAtom, direction, action, action);
    }

    @FunctionalInterface
    interface AtomPairConsumer extends BiConsumer<IAtom, IAtom> {

    }
}
