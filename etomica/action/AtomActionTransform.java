package etomica.action;
import etomica.*;
import etomica.space.Tensor;
import etomica.space.Vector;

public class AtomActionTransform extends AtomActionAdapter {
    
 //   private final AtomIteratorSequential iterator = new AtomIteratorSequential();
    
    public static void doAction(AtomIterator iterator, Vector r0, Tensor transformMatrix) {
        while(iterator.hasNext()) {
            iterator.nextAtom().coord.transform(r0, transformMatrix);
        }
    }
    
    public void actionPerformed(Atom a) {}
}