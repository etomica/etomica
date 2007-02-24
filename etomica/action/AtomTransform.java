package etomica.action;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.space.IVector;
import etomica.space.Tensor;

public class AtomTransform {
    
    public static void doTransform(AtomIterator iterator, IVector r0, Tensor transformMatrix) {
        while(iterator.hasNext()) {
            AtomLeaf a = (AtomLeaf)iterator.nextAtom();
            IVector r = a.getCoord().getPosition();
            r.ME(r0);
            a.getNode().parentPhase().getBoundary().nearestImage(r);
            r.transform(transformMatrix);
            r.PE(r0);
        }
    }
    
}