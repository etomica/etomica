package etomica;
import java.awt.*;
import java.awt.event.*;

/**
 * Superclass for classes that display information from simulation by painting to a canvas.
 * Defines methods useful for dealing with mouse and key events targeted at the display.
 * Much of the class is involved with defining event handling methods to permit display 
 * to be moved or resized; in the future these functions will be handled instead using awt component functions.
 * 
 * @see DisplayPhase.Canvas
 */
public abstract class DisplayCanvas extends javax.swing.JPanel implements java.io.Serializable, DisplayCanvasInterface, PhaseEventListener {

    protected Image offScreen;
    protected Graphics osg;
            
    protected DisplayPhase displayPhase;
    protected PhaseAction.Inflate inflate;

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
    
    /**
     *  Sets the quality of the rendered image, false = low, true = high
      */
    boolean highQuality = false;

    public DisplayCanvas() {
        setBackground(Color.white);
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
    
    public void phaseAction(PhaseEvent evt) {
        initialize();
    }
        
    public abstract void doPaint(Graphics g);
    
    public void update(Graphics g) {paint(g);}
        
    public void paint(Graphics g) {
        createOffScreen();
        doPaint(osg);
        g.drawImage(offScreen, 0, 0, null);
    }

    public void setPhase(Phase p) {
        inflate = new PhaseAction.Inflate(displayPhase.getPhase());
        p.speciesMaster.addListener(this);
    }
    
    public void setMovable(boolean b) {movable = b;}
    public boolean isMovable() {return movable;}
    public void setResizable(boolean b) {resizable = b;}
    public boolean isResizable() {return resizable;}

    public void setWriteScale(boolean s) {writeScale = s;}
    public boolean getWriteScale() {return(writeScale);}
    public void setHighQuality(boolean q) {highQuality = q;}
    public boolean getHighQuality() {return(highQuality);}

    public void initialize() {}

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

