package etomica.graphics;
import etomica.*;
import etomica.lattice.*;
import java.awt.*;

/**
 * Methods for drawing a lattice to a display.
 * Works for 2D only.
 */
 
 //plan to develop to configure for different types of drawing.
 //currently does only 2D cells
public class LatticeRenderer implements Drawable {
    
    public LatticeRenderer(AbstractLattice lattice) {
        this.lattice = lattice; 
        iterator.setBasis(lattice.siteList());
        
    }
    
    public void draw(java.awt.Graphics g, int[] origin, double scale) {
        g.setColor(Color.gray);
        double toPixels = scale*etomica.units.BaseUnit.Length.Sim.TO_PIXELS;
        iterator.reset();
        while(iterator.hasNext()) {
            AbstractCell cell = (AbstractCell)iterator.next();
            Space.Vector[] vertex = (Space.Vector[])cell.vertex();
            for(int i=1; i<vertex.length; i++) {
                Space2D.Vector v1 = null;
                Space2D.Vector v2 = null;
                switch(i) {//really specific!
                    case 0: 
                        v1 = (Space2D.Vector)vertex[0];
                        v2 = (Space2D.Vector)vertex[1];
                        break;
                    case 1: 
                        v1 = (Space2D.Vector)vertex[1];
                        v2 = (Space2D.Vector)vertex[3];
                        break;
                    case 2: 
                        v1 = (Space2D.Vector)vertex[2];
                        v2 = (Space2D.Vector)vertex[3];
                        break;
                    case 3: 
                        v1 = (Space2D.Vector)vertex[0];
                        v2 = (Space2D.Vector)vertex[2];
                        break;
                }
                int x1 = origin[0] + (int)(toPixels*v1.x);
                int y1 = origin[1] + (int)(toPixels*v1.y);
                int x2 = origin[0] + (int)(toPixels*v2.x);
                int y2 = origin[1] + (int)(toPixels*v2.y);
                g.drawLine(x1, y1, x2, y2);
            }
        }
        
    }
    
    private AbstractLattice lattice;
    private AtomIteratorListSimple iterator = new AtomIteratorListSimple();
}