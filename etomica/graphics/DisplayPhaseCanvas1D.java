package etomica.graphics;
import etomica.*;

import java.awt.*;
//import java.util.Hashtable;
import etomica.utility.java2.Iterator;

    /* History of changes
     * 7/16/02 (DAK) Modified for AtomType.Sphere diameter and radius method to take atom as argument.
     * 09/07/02 (DAK) added atomFilter
     */

//Class used to define canvas onto which configuration is drawn
public class DisplayPhaseCanvas1D extends DisplayCanvas {
    private TextField scaleText = new TextField();
    private Font font = new Font("sansserif", Font.PLAIN, 10);
    //  private int annotationHeight = font.getFontMetrics().getHeight();
    private int annotationHeight = 12;
    private int[] shiftOrigin = new int[2];     //work vector for drawing overflow images
    private final static Color wellColor = new Color(185,185,185, 110);
    private final AtomIteratorList atomIterator = new AtomIteratorList();
    private AtomFilter atomFilter;
    
    public DisplayPhaseCanvas1D(DisplayPhase _phase) {
        displayPhase = _phase;
        scaleText.setVisible(true);
        scaleText.setEditable(false);
        scaleText.setBounds(0,0,100,50);
        atomFilter = _phase.getAtomFilter();
    }
    
    public void initialize() {}
    
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
       
    private void drawAtom(Graphics g, int origin[], Atom a) {
        if(!atomFilter.accept(a)) return;
        
        Space.Vector r = a.coord.position();
        boolean drawWell = false;
        int sigmaP, xP, yP, baseXP, baseYP;

        if(a.type instanceof AtomTypeWell) {
            drawWell = true;
        }

        g.setColor(displayPhase.getColorScheme().atomColor(a));
            
        baseXP = origin[0] + (int)(displayPhase.getToPixels()*r.x(0));
        baseYP = origin[1] + Space1D.drawingHeight/2;
        if(a.type instanceof AtomTypeSphere) {
            /* Draw the core of the atom */
            sigmaP = (int)(displayPhase.getToPixels()*((AtomTypeSphere)a.type).diameter(a));
            xP = baseXP - (sigmaP>>1);
            yP = baseYP - (Space1D.drawingHeight >> 1);
            g.fillRect(xP, yP, sigmaP, Space1D.drawingHeight);
            /* Draw the surrounding well, if any */
            if(drawWell) {
                sigmaP = (int)(displayPhase.getToPixels()*((AtomTypeWell)a.type).wellDiameter());
                xP = baseXP - (sigmaP>>1);
                g.setColor(wellColor);
                g.drawRect(xP, yP, sigmaP, Space1D.drawingHeight);
            }
//            a.type.electroType().draw(g, origin, displayPhase.getToPixels(), r);
        } else if(a.type instanceof AtomTypeWall) {
            xP = baseXP;
            yP = baseYP - (Space1D.drawingHeight >> 1);
            int t = Math.max(1,(int)((double)((AtomTypeWall)a.type).getThickness()*(double)displayPhase.getToPixels()/(double)etomica.units.BaseUnit.Length.Sim.TO_PIXELS));
            if(!(((AtomTypeWall)a.type).isHorizontal() || ((AtomTypeWall)a.type).isVertical())) {  //not horizontal or vertical; draw line
                int x1 = xP + (int)(displayPhase.getToPixels()*((AtomTypeWall)a.type).getLength()*((AtomTypeWall)a.type).getCosX());
                int y1 = yP + (int)(displayPhase.getToPixels()*((AtomTypeWall)a.type).getLength()*((AtomTypeWall)a.type).getSinX());
                g.drawLine(xP, yP, x1, y1);
            }
            else if(((AtomTypeWall)a.type).isLongWall()) {
                int wP = ((AtomTypeWall)a.type).isVertical() ? t : Integer.MAX_VALUE;
                int hP = ((AtomTypeWall)a.type).isHorizontal() ? t : Space1D.drawingHeight;
                int X = ((AtomTypeWall)a.type).isVertical() ? xP : 0;
                int Y = ((AtomTypeWall)a.type).isHorizontal() ? yP : 0;
                g.fillRect(X,Y,wP,hP);
            }   
            else {                           //horizontal or vertical; draw box
                int wP = ((AtomTypeWall)a.type).isVertical() ? t : (int)(displayPhase.getToPixels()*((AtomTypeWall)a.type).getLength());
                int hP = ((AtomTypeWall)a.type).isHorizontal() ? t : Space1D.drawingHeight;
                g.fillRect(xP,yP,wP,hP);
            }
        } else { // Not a sphere, wall, or one of their derivatives...
            // Do nothing (how do you draw an object of unkown shape?)
        }
    }
            
    protected boolean computeShiftOrigin(Atom a, Space.Boundary b) {
        if(a.type instanceof AtomTypeSphere) {
            float[][] shifts = b.getOverflowShifts(a.coord.position(),((AtomTypeSphere)a.type).radius(a));  //should instead of radius have a size for all AtomC types
            for(int i=0; i<shifts.length; i++) {
                shiftOrigin[0] = displayPhase.getOrigin()[0] + (int)(displayPhase.getToPixels()*shifts[i][0]);
                shiftOrigin[1] = displayPhase.getOrigin()[1] + (int)(displayPhase.getToPixels()*shifts[i][1]);
            }
            return(true);
        } else {
            return(false);
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
            
        String vers = System.getProperty("java.version");
        if (vers.compareTo("1.2") >= 0) {
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
                RenderingHints.VALUE_ANTIALIAS_OFF);
            g2.setRenderingHint(RenderingHints.KEY_RENDERING, 
                RenderingHints.VALUE_RENDER_SPEED);
            g2.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, 
                RenderingHints.VALUE_COLOR_RENDER_SPEED);
            if(quality>=DRAW_QUALITY_HIGH) {
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
                    RenderingHints.VALUE_ANTIALIAS_ON);
                g2.setRenderingHint(RenderingHints.KEY_RENDERING, 
                    RenderingHints.VALUE_RENDER_QUALITY);
                g2.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, 
                    RenderingHints.VALUE_COLOR_RENDER_QUALITY);
            }
        }
            
        g.setColor(getBackground());
        g.fillRect(0,0,w,h);
        displayPhase.computeImageParameters2(w, h);

        //Draw other features if indicated
//        if(drawBoundary>DRAW_BOUNDARY_NONE) {displayPhase.getPhase().boundary().draw(g, displayPhase.getOrigin(), displayPhase.getScale());}
//        if(displayPhase.getDrawBoundary()) {displayPhase.getPhase().boundary().draw(g, displayPhase.getOrigin(), displayPhase.getScale());}

        //do drawing of all drawing objects that have been added to the display
        for(Iterator iter=displayPhase.getDrawables().iterator(); iter.hasNext(); ) {
            Drawable obj = (Drawable)iter.next();
            obj.draw(g, displayPhase.getOrigin(), displayPhase.getScale());
        }
            
        //Draw all atoms
        Space.Boundary boundary = displayPhase.getPhase().boundary();
        if(displayPhase.getColorScheme() instanceof ColorSchemeCollective) {
            ((ColorSchemeCollective)displayPhase.getColorScheme()).colorAllAtoms(displayPhase.getPhase());
        }
        atomIterator.setList(displayPhase.getPhase().speciesMaster.atomList);
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            drawAtom(g, displayPhase.getOrigin(), atomIterator.next());
        }
            
        //Draw overflow images if so indicated
        //This needs some work to make more general
///        for(Atom a=displayPhase.getPhase().firstAtom(); a!=null; a=a.nextAtom()) {
///            if(computeShiftOrigin(a, boundary)) {
///                drawAtom(g, shiftOrigin, a);
///            }
///        }

        //Draw periodic images if indicated
        if(displayPhase.getImageShells() > 0) {
            double[][] origins = displayPhase.getPhase().boundary().imageOrigins(displayPhase.getImageShells());  //more efficient to save rather than recompute each time
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
    }
}  //end of DisplayPhase.Canvas
