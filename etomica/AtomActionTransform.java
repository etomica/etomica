package etomica.action;
import etomica.*;

public class AtomActionTransform extends AtomAction {
    
    private final AtomIteratorSequential iterator = new AtomIteratorSequential();
    
    public static void doAction(AtomIterator iterator, Space.Vector r0, Space.Tensor rotationMatrix) {
        while(iterator.hasNext()) {
            iterator.next().coord.position().transform(r0, rotationMatrix);
        }
    }
    
    public void actionPerformed(Atom a) {}
}