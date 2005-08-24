package etomica.action;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIterator;
import etomica.space.Tensor;
import etomica.space.Vector;

public class AtomTransform {
    
    public static void doTransform(AtomIterator iterator, Vector r0, Tensor transformMatrix) {
        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            Vector r = a.coord.position();
            r.transform(a.node.parentPhase().boundary(),r0, transformMatrix);
        }
    }
    
}