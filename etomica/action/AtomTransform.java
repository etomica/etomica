package etomica.action;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.space.Tensor;
import etomica.space.Vector;

public class AtomTransform {
    
    public static void doTransform(AtomIterator iterator, Vector r0, Tensor transformMatrix) {
        while(iterator.hasNext()) {
            AtomLeaf a = (AtomLeaf)iterator.nextAtom();
            Vector r = a.coord.position();
            r.transform(a.getNode().parentPhase().getBoundary(),r0, transformMatrix);
        }
    }
    
}