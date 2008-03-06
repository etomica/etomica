package etomica.graphics;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.TextField;
import java.util.Iterator;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBoundary;
import etomica.api.IVector;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.atom.AtomTypeSphere;
import etomica.atom.AtomTypeWell;
import etomica.atom.IAtomOriented;
import etomica.space.Boundary;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Polygon;
import etomica.space.Space;

//Class used to define canvas onto which configuration is drawn
public class DisplayBoxCanvas2D extends DisplayCanvas {
    
    private TextField scaleText = new TextField();
    private Font font = new Font("sansserif", Font.PLAIN, 10);
    //  private int annotationHeight = font.getFontMetrics().getHeight();
    private int annotationHeight = 12;
    private int[] shiftOrigin = new int[2];     //work vector for drawing overflow images
    private final static Color wellColor = Color.pink;//new Color(185,185,185, 110);
    private final int[] atomOrigin;
    private final IVector boundingBox;
    protected final Space space;
        
    public DisplayBoxCanvas2D(DisplayBox _box, Space _space) {
    	this.space = _space;
        scaleText.setVisible(true);
        scaleText.setEditable(false);
        scaleText.setBounds(0,0,100,50);
        displayBox = _box;
        atomOrigin = new int[space.D()];
        boundingBox = space.makeVector();
    }
    
    /**
     * Sets the size of the display to a new value and scales the image so that
     * the box fits in the canvas in the same proportion as before.
     */
    public void scaleSetSize(int width, int height) {
        if(getBounds().width * getBounds().height != 0) {  //reset scale based on larger size change
            double ratio1 = (double)width/(double)getBounds().width;
            double ratio2 = (double)height/(double)getBounds().height;
            double factor = Math.min(ratio1, ratio2);
    //        double factor = (Math.abs(Math.log(ratio1)) > Math.abs(Math.log(ratio2))) ? ratio1 : ratio2;
            displayBox.setScale(displayBox.getScale()*factor);
            setSize(width, height);
        }
    }
          
    //Override superclass methods for changing size so that scale is reset with any size change  
    // this setBounds is ultimately called by all other setSize, setBounds methods
    public void setBounds(int x, int y, int width, int height) {
        if(width <= 0 || height <= 0) return;
        super.setBounds(x,y,width,height);
        createOffScreen(width,height);
    }
       
    protected void drawAtom(Graphics g, int origin[], IAtomPositioned a) {
        if(!displayBox.getAtomFilter().accept(a)) return;
        IVector r = a.getPosition();
        int sigmaP, xP, yP, baseXP, baseYP;

        boolean drawOrientation = (a.getType() instanceof AtomTypeOrientedSphere);
        boolean drawWell = (a.getType() instanceof AtomTypeWell);

        g.setColor(displayBox.getColorScheme().getAtomColor(a));

        baseXP = origin[0] + (int)(displayBox.getToPixels()*r.x(0));
        baseYP = origin[1] + (int)(displayBox.getToPixels()*r.x(1));
        if(a.getType() instanceof AtomTypeSphere) {
            /* Draw the core of the atom, specific to the dimension */
            sigmaP = (int)(displayBox.getToPixels()*((AtomTypeSphere)a.getType()).getDiameter());
            sigmaP = (sigmaP == 0) ? 1 : sigmaP;
            xP = baseXP - (sigmaP>>1);
            yP = baseYP - (sigmaP>>1);
            g.fillOval(xP, yP, sigmaP, sigmaP);
            /* Draw the surrounding well, if any, and specific to the dimension */
            if(drawWell) {
                sigmaP = (int)(displayBox.getToPixels()*((AtomTypeWell)a.getType()).wellDiameter());
                xP = baseXP - (sigmaP>>1);
                yP = baseYP - (sigmaP>>1);
                g.setColor(wellColor);
                g.drawOval(xP, yP, sigmaP, sigmaP);
            }
            /* Draw the orientation line, if any */
            if(drawOrientation) {
                IVector dir = ((IAtomOriented)a).getOrientation().getDirection();
                int dxy = (int)(displayBox.getToPixels()*0.5*((AtomTypeOrientedSphere)a.getType()).getDiameter());
                int dx = (int)(dxy*dir.x(0));
                int dy = (int)(dxy*dir.x(1));
                g.setColor(Color.red);
                xP += dxy; yP += dxy;
                g.drawLine(xP-dx, yP-dy, xP+dx, yP+dy);
            }
//            a.type.electroType().draw(g, origin, displayBox.getToPixels(), r);
        } else { // Not a sphere, wall, or one of their derivatives...
            // Do nothing (how do you draw an object of unknown shape?)
        }
    }
            
    IVector vec2;  
   /**
    * Method that handles the drawing of the box to the screen.
    *
    * @param g The graphic object to which the image of the box is drawn
    */
    public void doPaint(Graphics g) {
        if(!isVisible() || displayBox.getBox() == null) {return;}
        int w = getSize().width;
        int h = getSize().height;

        g.setColor(getBackground());
        g.fillRect(0,0,w,h);
        displayBox.computeImageParameters2(w, h);
        boundingBox.E(displayBox.getBox().getBoundary().getBoundingBox());
        int[] origin = displayBox.getOrigin();
        double toPixels = displayBox.getToPixels();

        //Draw other features if indicated
        // and the boundary is an etomica boundary.  Non-etomica objects
        // do not have a shape.
        if(drawBoundary>DRAW_BOUNDARY_NONE && displayBox.getBox().getBoundary() instanceof Boundary ) {
            g.setColor(Color.gray);
            Polygon shape = (Polygon)((Boundary)displayBox.getBox().getBoundary()).getShape();
            LineSegment[] edges = shape.getEdges();
            int ox = origin[0] + (int)(toPixels*boundingBox.x(0)*0.5);
            int oy = origin[1] + (int)(toPixels*boundingBox.x(1)*0.5);
            for(int i=0; i<edges.length; i++) {
                IVector[] vertices = edges[i].getVertices();
                int x1 = ox + (int)(toPixels*vertices[0].x(0));
                int y1 = oy + (int)(toPixels*vertices[0].x(1));
                int x2 = ox + (int)(toPixels*vertices[1].x(0));
                int y2 = oy + (int)(toPixels*vertices[1].x(1));
                g.drawLine(x1,y1,x2,y2);
            }
        }

//        if(displayBox.getDrawBoundary()) {displayBox.getBox().boundary().draw(g, displayBox.getOrigin(), displayBox.getScale());}

        //do drawing of all drawing objects that have been added to the display
        for(Iterator iter=displayBox.getDrawables().iterator(); iter.hasNext(); ) {
            Drawable obj = (Drawable)iter.next();
            obj.draw(g, origin, displayBox.getToPixels());
        }
            
        //Color all atoms according to colorScheme in DisplayBox
//        displayBox.getColorScheme().colorAllAtoms();

//        Vector[] vert = ((Polygon)displayBox.getBox().getBoundary().getShape()).getVertices();
//        vec2 = vert[2].M(vert[0]);
//        vec2.ME(vec);
//        vec2.TE(0.5);
        
        atomOrigin[0] = origin[0] + (int)(0.5*toPixels*boundingBox.x(0));
        atomOrigin[1] = origin[1] + (int)(0.5*toPixels*boundingBox.x(1));
//        vec.TE(0.0);
//        Vector vec2 = displayBox.getBox().getBoundary().centralImage(vec);

        //Draw all atoms
        if(displayBox.getColorScheme() instanceof ColorSchemeCollective) {
            ((ColorSchemeCollective)displayBox.getColorScheme()).colorAllAtoms();
        }
        IAtomSet leafList = displayBox.getBox().getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            if(this instanceof DisplayBoxSpin2D) {
            	drawAtom(g, origin, (IAtomPositioned)leafList.getAtom(iLeaf));
            }
            else {
                drawAtom(g, atomOrigin, (IAtomPositioned)leafList.getAtom(iLeaf));
            }
        }
            
        //Draw overflow images if so indicated
        if(displayBox.getDrawOverflow()) {
            IBoundary boundary = displayBox.getBox().getBoundary();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtomPositioned a = (IAtomPositioned)leafList.getAtom(iLeaf);
                if(!(a.getType() instanceof AtomTypeSphere)) continue;
                float[][] shifts = boundary.getOverflowShifts(a.getPosition(),0.5*((AtomTypeSphere)a.getType()).getDiameter());  //should instead of radius have a size for all AtomC types
                for(int i=shifts.length-1; i>=0; i--) {
                    shiftOrigin[0] = origin[0] + (int)(displayBox.getToPixels()*shifts[i][0]);
                    shiftOrigin[1] = origin[1] + (int)(displayBox.getToPixels()*shifts[i][1]);
                    drawAtom(g, shiftOrigin, a);
                }
            }
        }

        //Draw periodic images if indicated
        if(displayBox.getImageShells() > 0) {
            double[][] origins = displayBox.getBox().getBoundary().imageOrigins(displayBox.getImageShells());  //more efficient to save rather than recompute each time
            for(int i=0; i<origins.length; i++) {
                g.copyArea(displayBox.getOrigin()[0],displayBox.getOrigin()[1],displayBox.getDrawSize()[0],displayBox.getDrawSize()[1],(int)(displayBox.getToPixels()*origins[i][0]),(int)(displayBox.getToPixels()*origins[i][1]));
            }
        }
        //Draw bar showing scale if indicated
        if(writeScale) {
            g.setColor(Color.lightGray);
            g.fillRect(0,getSize().height-annotationHeight,getSize().width,annotationHeight);
            g.setColor(Color.black);
            g.setFont(font);
            g.drawString("Scale: "+Integer.toString((int)(100*displayBox.getScale()))+"%", 0, getSize().height-3);
        }
    }//end of doPaint
}  //end of DisplayBox.Canvas
