/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.ComponentListener;
import java.awt.event.KeyListener;
import java.awt.event.MouseListener;

import etomica.action.controller.Controller;
import etomica.units.Pixel;

/**
 * Superclass for classes that display information from simulation by painting to a canvas.
 * Defines methods useful for dealing with mouse and key events targeted at the display.
 * Much of the class is involved with defining event handling methods to permit display 
 * to be moved or resized; in the future these functions will be handled instead using awt component functions.
 * 
 */
public abstract class DisplayCanvas extends javax.swing.JPanel {

    public static final int DRAW_BOUNDARY_NONE = 0;
    public static final int DRAW_BOUNDARY_OUTLINE = 1;
    public static final int DRAW_BOUNDARY_SHELL = 2;
    public static final int DRAW_BOUNDARY_ALL = 3;
    public static final int DRAW_BOUNDARY_MAX = 4;

    protected Image offScreen;
    protected Graphics osg;
            
    protected DisplayBox displayBox;

    /**
    * Variable specifying whether a line tracing the boundary of the display should be drawn
    * Default value is <code>BOUNDARY_OUTLINE</code>
    */
    int drawBoundary = DRAW_BOUNDARY_OUTLINE;

    /**
     * Flag to indicate if display can be resized
     */
    boolean resizable = false;
    /**
     * Flag to indicate if display can be moved
     */
    boolean movable = false;

    /** 
     * Flag to indicate if value of scale should be superimposed on image
     */
    boolean writeScale = false;
    
    protected Pixel pixel;
    
    protected final Controller controller;
    
    /**
     * Construct a DisplayCanvas using the given controller (which may be null).
     * If a controller is given, the DisplayCanvas will suppress paint events
     * from the system that occur while the controller is active.
     */
    public DisplayCanvas(Controller controller) {
        this.controller = controller;
        setBackground(java.awt.Color.white);
    }

    public void dispose() {
        ComponentListener[] listeners = getComponentListeners();
        for (int i=0; i<listeners.length; i++) {
            removeComponentListener(listeners[i]);
        }
        MouseListener[] mlisteners = getMouseListeners();
        for (int i=0; i<mlisteners.length; i++) {
            removeMouseListener(mlisteners[i]);
        }
        KeyListener[] klisteners = getKeyListeners();
        for (int i=0; i<klisteners.length; i++) {
            removeKeyListener(klisteners[i]);
        }
    }

    protected void ensureOffScreen () {
        if (offScreen == null) { 
            createOffScreen(getSize().width, getSize().height);
        }
    }

    protected void createOffScreen(int w, int h) {
        offScreen = createImage(w,h);
        if(offScreen != null) osg = offScreen.getGraphics();
    }
    
    protected abstract void doPaint(Graphics g);
    
    public synchronized void paint(Graphics g) {
        if (!this.controller.isRunningActivityStep()) {
            // controller isn't running (we weren't called from the integrator)
            // so we need to do the drawing work here
            ensureOffScreen();
            if (osg == null) {
                return;
            }
            doPaint(osg);
        }
        // copy the image we just made (or was made previously on the
        // integrator's thread to the screen
        g.drawImage(offScreen, 0, 0, null);
    }

    public synchronized void repaint() {
        // do the drawing work now (on this thread)
        ensureOffScreen();
        if (osg == null) {
            return;
        }
        doPaint(osg);
        // now dispatch the paint request, which will happen on another thread
        super.repaint();
    }
    
    /**
     * Same as setSize, but included to implement DisplayCanvasInterface,
     * which has this for compatibility with OpenGL.
     */
    public void reshape(int width, int height) {
        setSize(width, height);
    }
    
    public void setMovable(boolean b) {movable = b;}
    public boolean isMovable() {return movable;}
    public void setResizable(boolean b) {resizable = b;}
    public boolean isResizable() {return resizable;}

    public void setWriteScale(boolean s) {writeScale = s;}
    public boolean getWriteScale() {return(writeScale);}

    public void setDrawBoundary(int b) {
      drawBoundary = b;
    }
    public int getDrawBoundary() {return drawBoundary;}

    public void setShiftX(float x) {}
    public void setShiftY(float y) {}
    public void setPrevX(float x) {}
    public void setPrevY(float y) {}
    public void setXRot(float x) {}
    public void setYRot(float y) {}
    public void setZoom(float z) {}
    public float getShiftX() {return(0f);}
    public float getShiftY() {return(0f);}
    public float getPrevX() {return(0f);}
    public float getPrevY() {return(0f);}
    public float getXRot() {return(0f);}
    public float getYRot() {return(0f);}
    public float getZoom() {return(1f);}
    public void startRotate(float x, float y) {};
    public void stopRotate() {};

    /**
     * Returns unit for conversion between simulation units and display pixels.
     */
    public Pixel getPixelUnit() {
        return pixel;
    }
    
    /**
     * Sets unit for conversion between simulation units and display pixels.
     */
    public void setPixelUnit(Pixel pixel) {
        this.pixel = pixel;
    }

} //end of DisplayCanvas class

