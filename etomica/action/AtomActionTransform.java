package etomica.action;
import etomica.Atom;
import etomica.AtomIterator;
import etomica.space.Boundary;
import etomica.space.Tensor;
import etomica.space.Vector;

public class AtomActionTransform extends AtomActionAdapter {
    
 //   private final AtomIteratorSequential iterator = new AtomIteratorSequential();
    
    public static void doAction(AtomIterator iterator, Vector r0, Tensor transformMatrix) {
        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            Vector r = a.coord.position();
            r.transform((Boundary)a.node.parentPhase().boundary(),r0, transformMatrix);
        }
    }
    
    public void actionPerformed(Atom a) {}
}