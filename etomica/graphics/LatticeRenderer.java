package etomica.graphics;
import etomica.*;
import etomica.lattice.*;
import java.awt.*;

/**
 * Methods for drawing a lattice to a display (2D only).
 * Defines color schemes useful to examination of cell neighbor list behavior.
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
                int x1 = origin[0] + (int)(toPixels*v1.x(0));
                int y1 = origin[1] + (int)(toPixels*v1.x(1));
                int x2 = origin[0] + (int)(toPixels*v2.x(0));
                int y2 = origin[1] + (int)(toPixels*v2.x(1));
                g.drawLine(x1, y1, x2, y2);
            }
        }
        
    }
    
    private AbstractLattice lattice;
    private AtomIteratorListSimple iterator = new AtomIteratorListSimple();
    
/////////////////////////////////////////////////////////////////////////////////////    
    
    public static class ColorSchemeNeighbor extends ColorSchemeCollective {
        
        private Atom referenceAtom;
        private final AtomIterator nbrIterator;
        private final AtomIteratorListSimple allIterator = new AtomIteratorListSimple();
        private final IteratorDirective directive = new IteratorDirective(IteratorDirective.BOTH);
    //    private final IteratorDirective directive = new IteratorDirective(IteratorDirective.UP);
        private final ColorSchemeByType typeColorScheme = new ColorSchemeByType();
        
        public ColorSchemeNeighbor(Simulation sim) {
            nbrIterator = sim.iteratorFactory.makeIntragroupNbrIterator();
        }
        
        public void colorAllAtoms(Phase phase) {
            allIterator.setBasis(phase.speciesMaster.atomList);
            allIterator.reset();
            while(allIterator.hasNext()) {
                Atom atom = allIterator.next();
                atom.allatomAgents[agentIndex] = typeColorScheme.atomColor(atom);//Color.green;
            }
            nbrIterator.reset(directive);
            while(nbrIterator.hasNext()) nbrIterator.next().allatomAgents[agentIndex] = Color.blue;
            referenceAtom.allatomAgents[agentIndex] = Color.red;
        }
        
        public void setAtom(Atom a) {
            referenceAtom = a;
            directive.set(a);
            nbrIterator.setBasis(a.node.parentGroup());
        }
        public Atom getAtom() {return referenceAtom;}
    }//end of ColorSchemeNeighbor
    
/////////////////////////////////////////////////////////////////////////////////////  

    public static class ColorSchemeCell extends ColorSchemeCollective {
        
        private int cellColorIndex = Atom.requestAgentIndex(this);
        private final AtomIteratorListSimple allIterator = new AtomIteratorListSimple();
        
        public void setLattice(AbstractLattice lattice) {
            AtomIteratorListSimple iterator = new AtomIteratorListSimple(lattice.siteList());
            while(iterator.hasNext()) {
                iterator.next().allatomAgents[cellColorIndex] = ConstantsGraphic.randomColor();
            }
        }
        
        public void colorAllAtoms(Phase phase) {
            allIterator.setBasis(phase.speciesMaster.atomList);
            allIterator.reset();
            while(allIterator.hasNext()) {
                Atom atom = allIterator.next();
                Atom cell = ((IteratorFactoryCell.NeighborSequencer)atom.seq).cell;
                atom.allatomAgents[agentIndex] = cell.allatomAgents[cellColorIndex];
            }
        }
    }//end of ColorSchemeCell
}//end of LatticeRenderer