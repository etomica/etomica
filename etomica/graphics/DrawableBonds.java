package etomica.graphics;
import etomica.*;
import java.awt.Graphics;

/**
 * Causes a line to be drawn between atoms joined by a Bond instance.
 * Appropriate only for 2D display.
 *
 * @author David Kofke
 */
 
 /* History
  * 12/04/02 (DAK) new
  */

public class DrawableBonds implements Drawable {
    
    private AtomIteratorListSimple iterator1;
    private AtomIteratorBonds iteratorBonds = new AtomIteratorBonds();
    private DisplayPhase displayPhase;
    private IteratorDirective directive = new IteratorDirective(IteratorDirective.DOWN);
    
    public DrawableBonds(DisplayPhase displayPhase, AtomList atoms) {
        iterator1 = new AtomIteratorListSimple(atoms);
        this.displayPhase = displayPhase;
    }
    
    public void draw(Graphics g, int[] origin, double s) {
        iterator1.reset();
        while(iterator1.hasNext()) {
            Atom atom1 = iterator1.next();
            Space.Vector r = atom1.coord.position();
            int x1 = origin[0] + (int)(displayPhase.getToPixels()*r.x(0));
            int y1 = origin[1] + (int)(displayPhase.getToPixels()*r.x(1));
      //      iteratorBonds.setBasis(atom1);
            iteratorBonds.reset(directive.set(atom1));
            while(iteratorBonds.hasNext()) {
                Atom atom2 = iteratorBonds.next();
                r = atom2.coord.position();
                int x2 = origin[0] + (int)(displayPhase.getToPixels()*r.x(0));
                int y2 = origin[1] + (int)(displayPhase.getToPixels()*r.x(1));
                g.drawLine(x1, y1, x2, y2);
            }
        }
    }
}
        
        