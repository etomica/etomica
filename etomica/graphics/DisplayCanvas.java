package etomica.graphics;
import java.awt.Graphics;
import java.awt.Image;

/**
 * Superclass for classes that display information from simulation by painting to a canvas.
 * Defines methods useful for dealing with mouse and key events targeted at the display.
 * Much of the class is involved with defining event handling methods to permit display 
 * to be moved or resized; in the future these functions will be handled instead using awt component functions.
 * 
 * @see DisplayPhase.Canvas
 */
public abstract class DisplayCanvas extends javax.swing.JPanel implements java.io.Serializable, DisplayCanvasInterface {

    protected Image offScreen;
    protected Graphics osg;
            
    protected DisplayPhase displayPhase;

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
     * Variable that sets the quality of the rendered image.
     */
    int quality = DRAW_QUALITY_NORMAL;

    /** 
     * Flag to indicate if value of scale should be superimposed on image
     */
    boolean writeScale = false;
    
    /**
     *  Sets the quality of the rendered image, false = low, true = high
      */
    boolean highQuality = false;

    public DisplayCanvas() {
        setBackground(java.awt.Color.white);
    }
    public void createOffScreen () {
        if (offScreen == null) { 
            createOffScreen(getSize().width, getSize().height);
        }
    }
    public void createOffScreen (int p) {
        createOffScreen(p,p);
    }
    public void createOffScreen(int w, int h) {
        offScreen = createImage(w,h);
        if(offScreen != null) osg = offScreen.getGraphics();
    }
    
    public abstract void doPaint(Graphics g);
    
    public void update(Graphics g) {paint(g);}
        
    public void paint(Graphics g) {
        createOffScreen();
        doPaint(osg);
        g.drawImage(offScreen, 0, 0, null);
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

    public void setQuality(int q) {
      if(q > DRAW_QUALITY_VERY_HIGH)
        q -= DRAW_QUALITY_MAX;
      if(q < DRAW_QUALITY_VERY_LOW)
        q += DRAW_QUALITY_MAX;
      quality = q;
    }
    public int getQuality() {return(quality);}
    
    public void setDrawBoundary(int b) {
      if(b>DRAW_BOUNDARY_ALL)
        b-=DRAW_BOUNDARY_MAX;
      else if(b<DRAW_BOUNDARY_NONE)
        b+=DRAW_BOUNDARY_MAX;
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

} //end of DisplayCanvas class

