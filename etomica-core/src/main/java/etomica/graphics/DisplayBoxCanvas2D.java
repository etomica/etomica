/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.action.controller.Controller;
import etomica.atom.*;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Polygon;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Pixel;

import java.awt.*;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.util.Iterator;

//Class used to define canvas onto which configuration is drawn
public class DisplayBoxCanvas2D extends DisplayCanvas {

    protected final Space space;
    private final int[] atomOrigin;
    private final Vector boundingBox;
    private Font font = new Font("sansserif", Font.PLAIN, 10);
    //  private int annotationHeight = font.getFontMetrics().getHeight();
    private int annotationHeight = 12;
    private int[] shiftOrigin = new int[2];     //work vector for drawing overflow images
    protected Color orientationColor = Color.RED;
    protected double orientationLength = 1.0;

    public DisplayBoxCanvas2D(DisplayBox _box, Space _space, Controller controller) {
        super(controller);
    	this.space = _space;

        displayBox = _box;
        atomOrigin = new int[space.D()];
        boundingBox = space.makeVector();

        addComponentListener(new ComponentListener() {
            public void componentHidden(ComponentEvent e) {}
            public void componentMoved(ComponentEvent e) {}
            public void componentShown(ComponentEvent e) {}
            public void componentResized(ComponentEvent e) { refreshSize(); }});
    }

    protected void refreshSize() {
        Dimension dim = getSize();
        Vector boxDim = displayBox.getBox().getBoundary().getBoxSize();
        double px = (dim.width - 1)/(boxDim.getX(0)+displayBox.getPaddingSigma());
        double py = (dim.height - 1)/(boxDim.getX(1)+displayBox.getPaddingSigma());
        if (px > py) {
            // take the smaller of the two toPixel
            px = py;
        }
        if (pixel != null && pixel.toPixels() == px) {
            return;
        }
        setPixelUnit(new Pixel(px));
        displayBox.computeImageParameters();
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

    protected void drawAtom(Graphics g, int[] origin, IAtom a) {
        Vector r = a.getPosition();
        int sigmaP, xP, yP, baseXP, baseYP;

        boolean drawOrientation = (a.getType() instanceof AtomTypeOriented);

        g.setColor(displayBox.getColorScheme().getAtomColor(a));

        double toPixels = pixel.toPixels() * displayBox.getScale();

        baseXP = origin[0] + (int)(toPixels*r.getX(0));
        baseYP = origin[1] + (int)(toPixels*r.getX(1));
        /* Draw the core of the atom, specific to the dimension */
        double sigma = displayBox.getDiameterHash().getDiameter(a);
        // default diameter
        if (sigma == -1) sigma = 1;
        sigmaP = (int)(toPixels*sigma);
        sigmaP = (sigmaP == 0) ? 1 : sigmaP;
        xP = baseXP - (sigmaP>>1);
        yP = baseYP - (sigmaP>>1);
        /* Draw the orientation line, if any */
        if(drawOrientation) {
            Vector dir = ((IAtomOriented)a).getOrientation().getDirection();
            double dxy = 2*(toPixels * 0.5 * sigma);
            int dx = (int) (.7*orientationLength * dxy * dir.getX(0));
            int dy = (int) (.7*orientationLength * dxy * dir.getX(1));
            g.setColor(Color.red);
            xP += (int) dxy;
            yP += (int) dxy;
            g.drawLine(xP, yP, xP + dx/2, yP + dy/2);
            g.fillOval(xP+dx/2-sigmaP/2, yP+dy/2-sigmaP/2, sigmaP, sigmaP);
            g.setColor(Color.black);
            g.drawLine(xP, yP, xP - dx/2, yP - dy/2);
            g.fillOval(xP-dx/2-sigmaP/2, yP-dy/2-sigmaP/2, sigmaP, sigmaP);
        } else {
            g.fillOval(xP, yP, sigmaP, sigmaP);
        }
    }

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
        boundingBox.E(displayBox.getBox().getBoundary().getBoxSize());
        int[] origin = displayBox.getOrigin();
        double toPixels = displayBox.getScale() * pixel.toPixels();

        //Draw other features if indicated
        // and the boundary is an etomica boundary.  Non-etomica objects
        // do not have a shape.
        if (drawBoundary > DRAW_BOUNDARY_NONE) {
            g.setColor(Color.gray);
            Polygon shape = (Polygon) displayBox.getBox().getBoundary().getShape();
            LineSegment[] edges = shape.getEdges();
            int ox = origin[0] + (int)(toPixels*boundingBox.getX(0)*0.5);
            int oy = origin[1] + (int)(toPixels*boundingBox.getX(1)*0.5);
            for(int i=0; i<edges.length; i++) {
                Vector[] vertices = edges[i].getVertices();
                int x1 = ox + (int)(toPixels*vertices[0].getX(0));
                int y1 = oy + (int)(toPixels*vertices[0].getX(1));
                int x2 = ox + (int)(toPixels*vertices[1].getX(0));
                int y2 = oy + (int)(toPixels*vertices[1].getX(1));
                g.drawLine(x1,y1,x2,y2);
            }
        }

//        if(displayBox.getDrawBoundary()) {displayBox.getBox().boundary().draw(g, displayBox.getOrigin(), displayBox.getScale());}


        //Color all atoms according to colorScheme in DisplayBox
//        displayBox.getColorScheme().colorAllAtoms();

//        Vector[] vert = ((Polygon)displayBox.getBox().getBoundary().getShape()).getVertices();
//        vec2 = vert[2].M(vert[0]);
//        vec2.ME(vec);
//        vec2.TE(0.5);

        atomOrigin[0] = origin[0] + (int)(0.5*toPixels*boundingBox.getX(0));
        atomOrigin[1] = origin[1] + (int)(0.5*toPixels*boundingBox.getX(1));
//        vec.TE(0.0);
//        Vector vec2 = displayBox.getBox().getBoundary().centralImage(vec);

        //Draw all atoms
        if(displayBox.getColorScheme() instanceof ColorSchemeCollective) {
            ((ColorSchemeCollective)displayBox.getColorScheme()).colorAllAtoms();
        }
        IAtomList leafList = displayBox.getBox().getLeafList();
        int nLeaf = leafList.size();
        AtomTest atomTest = displayBox.getAtomTestDoDisplay();
        if (atomTest instanceof AtomTestCollective) {
            ((AtomTestCollective) atomTest).resetTest();
        }
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom a = leafList.get(iLeaf);
            if (atomTest != null && !atomTest.test(a)) continue;
            drawAtom(g, atomOrigin, a);
        }

        //do drawing of all drawing objects that have been added to the display
        for(Iterator iter=displayBox.getDrawables().iterator(); iter.hasNext(); ) {
            Drawable obj = (Drawable)iter.next();
            obj.draw(g, origin, toPixels);
        }


        //Draw overflow images if so indicated
        if(displayBox.getDrawOverflow()) {
            Boundary boundary = displayBox.getBox().getBoundary();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtom a = leafList.get(iLeaf);
                OverflowShift overflow = new OverflowShift(space);
                double sigma = displayBox.getDiameterHash().getDiameter(a);
                if (sigma == -1) sigma = 1;
                float[][] shifts = overflow.getShifts(boundary, a.getPosition(),0.5*sigma);
                for(int i=shifts.length-1; i>=0; i--) {
                    shiftOrigin[0] = atomOrigin[0] + (int)(toPixels*shifts[i][0]);
                    shiftOrigin[1] = atomOrigin[1] + (int)(toPixels*shifts[i][1]);
                    drawAtom(g, shiftOrigin, a);
                }
            }
        }

        //Draw periodic images if indicated ONLY for an etomica Boundary
        if (displayBox.getImageShells() > 0) {

            double[][] origins = displayBox.getBox().getBoundary().imageOrigins(displayBox.getImageShells());  //more efficient to save rather than recompute each time
            for (int i = 0; i < origins.length; i++) {
                g.copyArea(displayBox.getOrigin()[0], displayBox.getOrigin()[1], displayBox.getDrawSize()[0], displayBox.getDrawSize()[1], (int) (toPixels * origins[i][0]), (int) (toPixels * origins[i][1]));
            }
        }
        //Draw bar showing scale if indicated
        if(writeScale) {
            g.setColor(Color.lightGray);
            g.fillRect(0,getSize().height-annotationHeight,getSize().width,annotationHeight);
            g.setColor(Color.black);
            g.setFont(font);
            g.drawString("Scale: " + (int) (100 * displayBox.getScale()) + "%", 0, getSize().height - 3);
        }
    }//end of doPaint

    public void setOrientationColor(Color color) {
        orientationColor = color;
    }

    /**
     * Sets the length of the drawn orientation line as a fraction of diameter
     *
     * @param orientationLength [copy/paste]
     */
    public void setOrientationLength(double orientationLength) {
        this.orientationLength = orientationLength;
    }
}  //end of DisplayBox.Canvas
