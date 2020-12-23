package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.NeighborConsumerHard;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManagerSimpleHard;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.collections.Int2IntHash;

/**
 * Simple neighbor iterator (iterates over all pairs) and can pass the pair
 * state to the callback for hard MD.
 */
public class NeighborIteratorSimpleHard implements NeighborIterator {
    private final NeighborManagerSimpleHard nbrManager;
    private final Box box;
    private final Space space;

    public NeighborIteratorSimpleHard(NeighborManagerSimpleHard nbrManager, Box box) {
        this.nbrManager = nbrManager;
        this.box = box;
        this.space = box.getSpace();
    }

    @Override
    public void iterUpNeighbors(int i, NeighborConsumer consumer) {
        final NeighborConsumerHard consumerHard = consumer instanceof NeighborConsumerHard ? (NeighborConsumerHard) consumer : null;

        IAtomList atoms = box.getLeafList();
        IAtom atom1 = atoms.get(i);
        Vector ri = atom1.getPosition();
        Int2IntHash iState = nbrManager.stateHash[i];

        for (int j = i + 1; j < atoms.size(); j++) {
            IAtom jAtom = atoms.get(j);
            Vector rj = jAtom.getPosition();
            Vector rij = space.makeVector();
            rij.Ev1Mv2(rj, ri);
            box.getBoundary().nearestImage(rij);
            if (consumerHard == null) {
                consumer.accept(jAtom, rij);
            } else {
                consumerHard.accept(jAtom, rij, iState.get(j));
            }
        }
    }

    @Override
    public void iterDownNeighbors(int i, NeighborConsumer consumer) {
        final NeighborConsumerHard consumerHard = consumer instanceof NeighborConsumerHard ? (NeighborConsumerHard) consumer : null;

        IAtomList atoms = box.getLeafList();
        IAtom atom1 = atoms.get(i);
        Vector ri = atom1.getPosition();

        for (int j = 0; j < i; j++) {
            IAtom jAtom = atoms.get(j);
            Vector rj = jAtom.getPosition();
            Vector rij = space.makeVector();
            rij.Ev1Mv2(rj, ri);
            box.getBoundary().nearestImage(rij);
            if (consumerHard == null) {
                consumer.accept(jAtom, rij);
            } else {
                consumerHard.accept(jAtom, rij, nbrManager.stateHash[j].get(i));
            }
        }
    }

    @Override
    public void iterAllNeighbors(int iAtom, NeighborConsumer consumer) {
        throw new UnsupportedOperationException();
    }
}
