package etomica.action;
import etomica.*;

public class AtomActionTransform extends AtomAction {
    
 //   private final AtomIteratorSequential iterator = new AtomIteratorSequential();
    
    public static void doAction(AtomIterator iterator, Space.Vector r0, Space.Tensor transformMatrix) {
        while(iterator.hasNext()) {
            iterator.nextAtom().coord.transform(r0, transformMatrix);
        }
    }
    
    public void actionPerformed(Atom a) {}
}