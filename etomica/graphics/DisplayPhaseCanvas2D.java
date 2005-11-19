package etomica.graphics;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.TextField;
import java.util.Iterator;

import etomica.atom.Atom;
import etomica.atom.AtomFilter;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.atom.AtomTypeSphere;
import etomica.atom.AtomTypeWell;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Polygon;
import etomica.space.Boundary;
import etomica.space.ICoordinateAngular;
import etomica.space.Vector;
import etomica.species.Species;

    /* History of changes
     * 7/16/02 (DAK) Modified for AtomType.Sphere diameter and radius method to take atom as argument.
     * 8/18/02 (DAK) Modified drawing of Wall to shift image by drawShift value given in type
     * 9/07/02 (DAK) Added atomFilter
     */

//Class used to define canvas onto which configuration is drawn
public class DisplayPhaseCanvas2D extends DisplayCanvas {
    
    private TextField scaleText = new TextField();
    private Font font = new Font("sansserif", Font.PLAIN, 10);
    //  private int annotationHeight = font.getFontMetrics().getHeight();
    private int annotationHeight = 12;
    private int[] shiftOrigin = new int[2];     //work vector for drawing overflow images
    private final static Color wellColor = Color.pink;//new Color(185,185,185, 110);
    private final AtomIteratorLeafAtoms atomIterator = new AtomIteratorLeafAtoms();
    private AtomFilter atomFilter;
        
    public DisplayPhaseCanvas2D(DisplayPhase _phase) {
        scaleText.setVisible(true);
        scaleText.setEditable(false);
        scaleText.setBounds(0,0,100,50);
        displayPhase = _phase;
        atomFilter = _phase.getAtomFilter();
    }
    
    public void setAtomFilter(AtomFilter filter) {atomFilter = filter;}
    
    /**
     * Sets the size of the display to a new value and scales the image so that
     * the phase fits in the canvas in the same proportion as before.
     */
    public void scaleSetSize(int width, int height) {
        if(getBounds().width * getBounds().height != 0) {  //reset scale based on larger size change
            double ratio1 = (double)width/(double)getBounds().width;
            double ratio2 = (double)height/(double)getBounds().height;
            double factor = Math.min(ratio1, ratio2);
    //        double factor = (Math.abs(Math.log(ratio1)) > Math.abs(Math.log(ratio2))) ? ratio1 : ratio2;
            displayPhase.setScale(displayPhase.getScale()*factor);
            setSize(width, height);
        }
    }
          
    //Override superclass methods for changing size so that scale is reset with any size change  
    // this setBounds is ultimately called by all other setSize, setBounds methods
    public void setBounds(int x, int y, int width, int height) {
        if(width == 0 || height == 0) return;
        super.setBounds(x,y,width,height);
        createOffScreen(width,height);
    }
       
    protected void drawAtom(Graphics g, int origin[], Atom a) {
        if(!atomFilter.accept(a)) return;
        Vector r = a.coord.position();
        int sigmaP, xP, yP, baseXP, baseYP;

        boolean drawOrientation = (a.type instanceof AtomTypeOrientedSphere);
        boolean drawWell = (a.type instanceof AtomTypeWell);

        g.setColor(displayPhase.getColorScheme().atomColor(a));
            
        baseXP = origin[0] + (int)(displayPhase.getToPixels()*r.x(0));
        baseYP = origin[1] + (int)(displayPhase.getToPixels()*r.x(1));
        if(a.type instanceof AtomTypeSphere) {
            /* Draw the core of the atom, specific to the dimension */
            sigmaP = (int)(displayPhase.getToPixels()*((AtomTypeSphere)a.type).diameter(a));
            sigmaP = (sigmaP == 0) ? 1 : sigmaP;
            xP = baseXP - (sigmaP>>1);
            yP = baseYP - (sigmaP>>1);
            g.fillOval(xP, yP, sigmaP, sigmaP);
            /* Draw the surrounding well, if any, and specific to the dimension */
            if(drawWell) {
                sigmaP = (int)(displayPhase.getToPixels()*((AtomTypeWell)a.type).wellDiameter());
                xP = baseXP - (sigmaP>>1);
                yP = baseYP - (sigmaP>>1);
                g.setColor(wellColor);
                g.drawOval(xP, yP, sigmaP, sigmaP);
            }
            /* Draw the orientation line, if any */
            if(drawOrientation) {
                double theta = ((ICoordinateAngular)a.coord).orientation().angle()[0];
                int dxy = (int)(displayPhase.getToPixels()*((AtomTypeOrientedSphere)a.type).radius(a));
                int dx = (int)(dxy*Math.cos(theta));
                int dy = (int)(dxy*Math.sin(theta));
                g.setColor(Color.red);
                xP += dxy; yP += dxy;
                g.drawLine(xP-dx, yP-dy, xP+dx, yP+dy);
            }
//            a.type.electroType().draw(g, origin, displayPhase.getToPixels(), r);
        } else { // Not a sphere, wall, or one of their derivatives...
            // Do nothing (how do you draw an object of unknown shape?)
        }
    }
            
      
    /**
    * doPaint is the method that handles the drawing of the phase to the screen.
    * Several variables and conditions affect how the image is drawn.  First,
    * the Unit.Length.Sim class variable <code>TO_PIXELS</code> performs the conversion between simulation
    * length units (Angstroms) and pixels.  The default value is 10 pixels/Angstrom
    * reflecting the default size of the phase (300 pixels by 300 pixels) and the
    * default phase size (30 by 30 A).  
    * The field <code>scale</code> is a multiplicative factor that directly
    * scales up or down the size of the image; this value is adjusted automatically
    * whenever shells of periodic images are drawn, to permit the central image and all 
    * of the specified periodic images to fit in the drawing of the phase.  
    *
    * @param g The graphic object to which the image of the phase is drawn
    * @see Species
    */
    public void doPaint(Graphics g) {
        if(!isVisible() || displayPhase.getPhase() == null) {return;}
        int w = getSize().width;
        int h = getSize().height;

        g.setColor(getBackground());
        g.fillRect(0,0,w,h);
        displayPhase.computeImageParameters2(w, h);

        //Draw other features if indicated
        if(drawBoundary>DRAW_BOUNDARY_NONE) {
            g.setColor(Color.gray);
            double toPixels = displayPhase.getScale()*etomica.units.BaseUnit.Length.Sim.TO_PIXELS;
            Polygon shape = (Polygon)displayPhase.getPhase().getBoundary().getShape();
            LineSegment[] edges = shape.getEdges();
            int ox = displayPhase.getOrigin()[0] + (int)(toPixels*displayPhase.getPhase().getBoundary().getDimensions().x(0)*0.5);
            int oy = displayPhase.getOrigin()[1] + (int)(toPixels*displayPhase.getPhase().getBoundary().getDimensions().x(1)*0.5);
            for(int i=0; i<edges.length; i++) {
                int x1 = ox + (int)(toPixels*edges[i].getVertices()[0].x(0));
                int y1 = oy + (int)(toPixels*edges[i].getVertices()[0].x(1));
                int x2 = ox + (int)(toPixels*edges[i].getVertices()[1].x(0));
                int y2 = oy + (int)(toPixels*edges[i].getVertices()[1].x(1));
                g.drawLine(x1,y1,x2,y2);
            }
        }


//        if(displayPhase.getDrawBoundary()) {displayPhase.getPhase().boundary().draw(g, displayPhase.getOrigin(), displayPhase.getScale());}

        //do drawing of all drawing objects that have been added to the display
        for(Iterator iter=displayPhase.getDrawables().iterator(); iter.hasNext(); ) {
            Drawable obj = (Drawable)iter.next();
            obj.draw(g, displayPhase.getOrigin(), displayPhase.getScale());
        }
            
        //Color all atoms according to colorScheme in DisplayPhase
//        displayPhase.getColorScheme().colorAllAtoms();
            
        //Draw all atoms
        Boundary boundary = displayPhase.getPhase().getBoundary();
        if(displayPhase.getColorScheme() instanceof ColorSchemeCollective) {
            ((ColorSchemeCollective)displayPhase.getColorScheme()).colorAllAtoms(displayPhase.getPhase());
        }
        atomIterator.setPhase(displayPhase.getPhase());
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            drawAtom(g, displayPhase.getOrigin(), atomIterator.nextAtom());
        }
            
        //Draw overflow images if so indicated
        if(displayPhase.getDrawOverflow()) {
            atomIterator.reset();
            while(atomIterator.hasNext()) {
                Atom a = atomIterator.nextAtom();
                if(!(a.type instanceof AtomTypeSphere)) continue;
                float[][] shifts = boundary.getOverflowShifts(a.coord.position(),((AtomTypeSphere)a.type).radius(a));  //should instead of radius have a size for all AtomC types
                for(int i=shifts.length-1; i>=0; i--) {
                    shiftOrigin[0] = displayPhase.getOrigin()[0] + (int)(displayPhase.getToPixels()*shifts[i][0]);
                    shiftOrigin[1] = displayPhase.getOrigin()[1] + (int)(displayPhase.getToPixels()*shifts[i][1]);
                    drawAtom(g, shiftOrigin, a);
                }
            }
        }

        //Draw periodic images if indicated
        if(displayPhase.getImageShells() > 0) {
            double[][] origins = displayPhase.getPhase().getBoundary().imageOrigins(displayPhase.getImageShells());  //more efficient to save rather than recompute each time
            for(int i=0; i<origins.length; i++) {
                g.copyArea(displayPhase.getOrigin()[0],displayPhase.getOrigin()[1],displayPhase.getDrawSize()[0],displayPhase.getDrawSize()[1],(int)(displayPhase.getToPixels()*origins[i][0]),(int)(displayPhase.getToPixels()*origins[i][1]));
            }
        }
        //Draw bar showing scale if indicated
        if(writeScale) {
            g.setColor(Color.lightGray);
            g.fillRect(0,getSize().height-annotationHeight,getSize().width,annotationHeight);
            g.setColor(Color.black);
            g.setFont(font);
            g.drawString("Scale: "+Integer.toString((int)(100*displayPhase.getScale()))+"%", 0, getSize().height-3);
        }
    }//end of doPaint
}  //end of DisplayPhase.Canvas
