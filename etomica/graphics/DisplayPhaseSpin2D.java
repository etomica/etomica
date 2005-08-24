package etomica.graphics;

import java.awt.Color;
import java.awt.Graphics;

import etomica.atom.Atom;
import etomica.lattice.RectangularLattice;
import etomica.nbr.site.AtomSequencerSite;
import etomica.nbr.site.AtomSite;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 23, 2005 by kofke
 */
public class DisplayPhaseSpin2D extends DisplayPhaseCanvas2D {

    /**
     * @param phase
     */
    public DisplayPhaseSpin2D(DisplayPhase _phase) {
        super(_phase);
        latticeIndex = new int[displayPhase.getPhase().space().D()];
        lattice = displayPhase.getPhase().getLattice();
        spinWidth = 5;
    }
    
    protected void drawAtom(Graphics g, int origin[], Atom atom) {
        AtomSite site = ((AtomSequencerSite)atom.seq).getSite();
        if (site == null) return;
        displayPhase.getPhase().getLattice().latticeIndex(site.getLatticeArrayIndex(),latticeIndex);
            
        //color central site red
//        ((MySite)lattice.site(iterator.centralSite)).color = Color.RED;
        int nx = 5; //number of planes to draw before moving down to another line of planes
        int k = latticeIndex.length < 3 ? 0 : latticeIndex[2];
        int ox = origin[0] + k % nx;       //set origin for drawing each lattice plane
        int oy = origin[1] + (k - ox)/nx;
        ox = 3 + ox*spinWidth*(lattice.getSize()[0]+2);
        oy = 3 + oy*spinWidth*(lattice.getSize()[1]+2);
        //draw lattice plane
        g.setColor(atom.coord.position().x(0) > 0 ? Color.green : Color.white);
        g.fillRect(ox+latticeIndex[0]*spinWidth,oy+latticeIndex[1]*spinWidth,spinWidth,spinWidth);
//        g.setColor(Color.black);
//        g.drawRect(ox+latticeIndex[0]*spinWidth,oy+latticeIndex[1]*spinWidth,spinWidth,spinWidth);
    }

    private int spinWidth;
    private final int[] latticeIndex;
    private final RectangularLattice lattice;
}
