package etomica;
import java.awt.*;
import java.awt.event.*;
import java.util.Iterator; 
import java.util.Observer;
import java.util.Observable;
import java.util.Vector;
import etomica.units.*;
import javax.swing.JDialog;

/**
 * Displays a picture of a phase, with configurations of molecules, boundaries, and other objects as appropriate, assuming 2-dimensional system.  
 * Instantiates a ConfigurationCanvas (an inner class of this one) for most of the work.
 * DisplayPhase is an input event (mouse and key) listener of the canvas.  It receives these 
 * events and uses information from them to form and fire a DisplayPhaseEvent to registered listeners.
 *
 * @see DisplayPhaseEvent
 * @author David Kofke
 */
public class DisplayPhase extends Display {
        
    public static final int LEFT = -1;   //Class variables to code for alignment of drawn image within display region
    public static final int CENTER = 0;
    public static final int RIGHT = +1;
    public static final int TOP = -1;
    public static final int BOTTOM = +1;
    private final int D = 2;
    protected ColorScheme colorScheme = new ColorSchemeByType();
        
//    private final ConfigurationCanvas canvas = new ConfigurationCanvas();
    public Canvas canvas;  //do not instantiate here; instead must be in graphic method

//Explicit to 2D because drawing to 2D image
    public final int[] align = new int[D];
    
 /**
  * Size of drawing region of central image, in pixels
  *
  * @see #computeDrawSize
  */
    protected final int[] drawSize = new int[D];
  
 /**
  * Flag specifying whether a line tracing the boundary of the display should be drawn
  * Default value is <code>true</code>
  */
  private boolean drawBoundary = true;
  
 /**
  * Flag specifying whether a line tracing the boundary of the space should be drawn
  * Default value is <code>false</code>
  */
  private boolean drawSpace = false;
  
 /**
  * Flag specifying whether meters should be allowed to draw something
  * Default value is <code>false</code>
  */
  private boolean drawMeters = false;

 /**
  * Number of periodic-image shells to be drawn when drawing this phase to the
  * screen.  Default value is 0.
  *
  * @see #paint
  */
  private int imageShells = 0;
 
 /**
  * Factor used to scale the size of the image. May be used
  * to scale up or down the image within one phase without affecting those
  * in other displays.  Default value is 1.0.
  */
    protected double scale = 1.0;
      
   /**
    * Coordinate origin for central image
    * Explicit to 2D because drawing is done to 2D image
    */
    protected final int[] centralOrigin = new int[D];

    private transient final int[] shiftOrigin = new int[D];     //work vector for drawing overflow images

 /**
  * When using periodic boundaries, image molecules near the cell boundaries often have parts that overflow
  * into the central cell.  When the phase is drawn, these "overflow portions" are not normally
  * included in the central image.  Setting this flag to <code>true</code> causes extra drawing
  * to be done so that the overflow portions are properly rendered.  This is particularly helpful
  * to have on when imageShells is non-zero.  Default value is <code>false</code>.
  */
  private boolean drawOverflow = false;
  
  /**
   * Vector used to maintain list of DisplayPhase Listeners
   */
  private java.util.Vector displayPhaseListeners = new java.util.Vector();
  
  /**
   * Iterator of atoms in the displayed phase
   */
   private Atom.Iterator atomIterator;
  
    public DisplayPhase () {
        this(Simulation.instance);
    }
    public DisplayPhase(Simulation sim) {
        super(sim);
        
        setLabel("Configuration");

        canvas = new DisplayPhase.Canvas();

        int box = (int)(Default.BOX_SIZE * BaseUnit.Length.Sim.TO_PIXELS);
        setSize(box, box);
        align[0] = align[1] = CENTER;
        
        InputEventHandler listener = new InputEventHandler();
        canvas.addMouseListener((MouseListener)listener);
        canvas.addMouseMotionListener((MouseMotionListener)listener);
        canvas.addKeyListener((KeyListener)listener);
        
        canvas.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent evt) {
                if((evt.getModifiers() & InputEvent.BUTTON3_MASK) != 0) {
                    if(DeviceConfigurationEditor.exists) return;
                    Device editor = new DeviceConfigurationEditor(DisplayPhase.this);
                    parentSimulation().add(editor.graphic(null));
                    parentSimulation().validate();
                    parentSimulation().repaint();
                }
            }
        });

        setLayout(null);
    }
    
    public void setBounds(int x, int y, int width, int height) {
        super.setBounds(x,y,width,height);
        canvas.setBounds(x,y,width,height);
    }

    public void setAlign(int i, int value) {
        align[i] = value;
    }
    public int getAlign(int i) {return align[i];}

    public final boolean getDrawOverflow() {return drawOverflow;}
    public final void setDrawOverflow(boolean b) {drawOverflow = b;}

    public double getScale() {return scale;}
    public void setScale(double s) {
        if(s>0) {
            scale = s;
        }
    }
    
    /**
     * Specifies the phase for this display
     */
    public void setPhase(Phase p) {
        if(p == null) return;
        super.setPhase(p);
        colorScheme.setPhase(p);
        canvas.setPhase(p);
    }//need to add an iterator observer (in superclass)
    
    /**
     * Accessor method for the color scheme used for this display
     */
    public void setColorScheme(ColorScheme colorScheme) {
        this.colorScheme = colorScheme;
        if(phase != null) colorScheme.setPhase(phase);
    }
    /**
     * Accessor method for the color scheme used for this display
     */
    public ColorScheme getColorScheme() {return colorScheme;}
    
    /** 
     * Simulation.GraphicalElement interface method.  Overrides Display method
     * to return the DisplayPhase.Canvas as the display object.
     *
     * @param obj ignored by this method.
     */
    public Component graphic(Object obj) {
        return canvas;
    }

    /**
    * @return the current value of imageShells.
    */
    public int getImageShells() {return imageShells;}
     
    /**
    * Changes the value of image shells, and increases/decreases scale accordingly.
    *
    * @param n the new value of imageShells
    */
    public void setImageShells(int n) {
        if(n>=0) {
            scale *= (double)(2*imageShells+1)/(double)(2*n+1);
            imageShells = n;
        }
    }
      
    protected int computeOrigin(int align, int drawSize, int size) {
        switch(align) {
            case   LEFT: return 0;    //same as TOP
            case CENTER: return (size-drawSize)/2;
            case  RIGHT: return size-drawSize; //same as BOTTOM
            default: return 0;
        }
    }
            
    public void setDrawBoundary(boolean b) {drawBoundary = b;}
    public boolean getDrawBoundary() {return drawBoundary;}
    public void setDrawSpace(boolean b) {drawSpace = b;}
    public boolean getDrawSpace() {return drawSpace;}
    public void setDrawMeters(boolean b) {drawMeters = b;}
    public boolean getDrawMeters() {return drawMeters;}

    public void doUpdate() {;}
    public void repaint() {canvas.repaint();}
      
    public void setMovable(boolean b) {canvas.setMovable(b);}
    public boolean isMovable() {return canvas.isMovable();}
    public void setResizable(boolean b) {canvas.setResizable(b);}
    public boolean isResizable() {return canvas.isResizable();}
    
    //Methods for handling DisplayPhaseEvents
    
    public synchronized void addDisplayPhaseListener(DisplayPhaseListener dpl) {
        displayPhaseListeners.addElement(dpl);
    }
    public synchronized void removeDisplayPhaseListener(DisplayPhaseListener dpl) {
        displayPhaseListeners.removeElement(dpl);
    }
    public void fireDisplayPhaseEvent(DisplayPhaseEvent dpe) {
//        Vector currentListeners = null;
//        synchronized(this){
//            currentListeners = (Vector)displayPhaseEventListeners.clone();
//        }
//        for(int i = 0; i < currentListeners.size(); i++) {
//            DisplayPhaseEventListener listener = (DisplayPhaseEventListener)currentListeners.elementAt(i);
//            listener.displayPhaseEventAction(dpe);
//        }
        for(int i = 0; i < displayPhaseListeners.size(); i++) {
            DisplayPhaseListener listener = (DisplayPhaseListener)displayPhaseListeners.elementAt(i);
            listener.displayPhaseAction(dpe);
        }
    }
    

    //Class used to define canvas onto which configuration is drawn
    
    public class Canvas extends DisplayCanvas {
        
        private boolean writeScale = false;  //flag to indicate if value of scale should be superimposed on image
        private TextField scaleText = new TextField();
        private Font font = new Font("sansserif", Font.PLAIN, 10);
        //  private int annotationHeight = font.getFontMetrics().getHeight();
        private int annotationHeight = 12;
        protected PhaseAction.Inflate inflate;
        
        public Canvas() {
            scaleText.setVisible(true);
            scaleText.setEditable(false);
            scaleText.setBounds(0,0,100,50);
        }
        
        public void setPhase(Phase p) {
            inflate = new PhaseAction.Inflate(phase);
        }
        
        public void changeSize(int w, int h, MouseEvent e) {
            if(e.isControlDown()) {  //change volume of simulated phase
                double rScale;
                int wh;
                if(w > h) {
                    wh = w;
                    rScale = (double)w/(double)getSize().width;
                }
                else {
                    wh = h;
                    rScale = (double)h/(double)getSize().height;
                }
                super.setSize(wh,wh);
            //    createOffScreen(wh);  //redundant with call now in setBounds
                inflate.actionPerformed(phase,rScale);
                phase.integrator().initialize();
            }
            else {                    //change scale of image
                //super.changeSize(w,h,e);  doesn't work well
                int x = Math.max(w,h);
                setSize(x,x);
             //   createOffScreen(x);  //redundant
            }
        }
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
                scale *= factor;
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
       
        /**
         * Defines action of display so that it calls setImageCells when a digit is typed.
         */
        public void digitTyped(int i) {
            setImageShells(i);
        }
        /**
         * Performs actions when some particular letters are typed.  In particular<ul>
         * <li>"s" toggles presentation of image scale to screen
         * <li>"o" toggles drawing of "overflow images" (parts of molecules from surrounding periodic images
         * <li>"b" toggles drawing of boundary
         * </ul>
         */
        public void letterTyped(char c) {
            switch(c) {
                case 's':
                    writeScale = !writeScale;
                    break;
                case 'o':
                    drawOverflow = !drawOverflow;
                    break;
                case 'b':
                    setDrawBoundary(!getDrawBoundary());
                    break;
                default:
                    break;
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
        * Main method for performing the paint is in each atom.
        *
        * @param g The graphic object to which the image of the phase is drawn
        * @see Species
        */
        public void doPaint(Graphics g) {  //1-D or 2-D
            if(!isVisible() || phase() == null) {return;}
            int w = getSize().width;
            int h = getSize().height;
            g.setColor(getBackground());
            g.fillRect(0,0,w,h);
            //Compute factor converting simulation units to pixels for this display
            double toPixels = getScale()*BaseUnit.Length.Sim.TO_PIXELS;
            //Determine length and width of drawn image, in pixels
            drawSize[0] = (int)(toPixels*phase().boundary().dimensions().component(0));
            drawSize[1] = (parentSimulation().space().D()==1) ? Space1D.drawingHeight: (int)(toPixels*phase().boundary().dimensions().component(1));
            //Find origin for drawing action
            centralOrigin[0] = computeOrigin(align[0],drawSize[0],w);
            centralOrigin[1] = computeOrigin(align[1],drawSize[1],h);
            //Draw other features if indicated
            if(drawBoundary) {phase().boundary().draw(g, centralOrigin, scale);}
            if(drawSpace) {parentSimulation().space().draw(g, centralOrigin, scale);}
//            if(drawMeters) {
//                for(java.util.Iterator iter=phase.meterManager().list.iterator(); iter.hasNext(); ) {
//                    ((MeterAbstract)iter.next()).draw(g, centralOrigin, scale);
//                }
//            }
            //Color all atoms according to colorScheme in DisplayPhase
            colorScheme.colorAllAtoms();
            
            //Draw all atoms
            Space.Boundary boundary = phase().boundary();
            for(Atom a = phase().firstAtom(); a!=null; a=a.nextAtom()) {
                boundary.centralImage(a.coordinate);        //move atom to central image
                a.draw(g,centralOrigin,toPixels);
            }
            //Draw overflow images if so indicated
            //This needs some work to make more general
            if(drawOverflow) {
                for(Atom a=phase().firstAtom(); a!=null; a=a.nextAtom()) {
                    if(a.type instanceof AtomType.Disk) {
                        double[][] shifts = boundary.getOverflowShifts(a.coordinate.position(),((AtomType.Disk)a.type).radius());  //should instead of radius have a size for all AtomC types
                        for(int i=0; i<shifts.length; i++) {
                        shiftOrigin[0] = centralOrigin[0] + (int)(toPixels*shifts[i][0]);
                        shiftOrigin[1] = centralOrigin[1] + (int)(toPixels*shifts[i][1]);
                        a.draw(g,shiftOrigin,toPixels);
                        }
                    }
                }
            } 
            //Draw periodic images if indicated
            if(imageShells > 0) {
                double[][] origins = phase().boundary().imageOrigins(imageShells);  //more efficient to save rather than recompute each time
                for(int i=0; i<origins.length; i++) {
                    g.copyArea(centralOrigin[0],centralOrigin[1],drawSize[0],drawSize[1],(int)(toPixels*origins[i][0]),(int)(toPixels*origins[i][1]));
                }
            }
            //Draw bar showing scale if indicated
            if(writeScale) {
                g.setColor(Color.lightGray);
                g.fillRect(0,getSize().height-annotationHeight,getSize().width,annotationHeight);
                g.setColor(Color.black);
                g.setFont(font);
                g.drawString("Scale: "+Integer.toString((int)(100*scale))+"%", 0, getSize().height-3);
            }
        }
    }  //end of DisplayConfiguration.Canvas
    
    /**
     * Class to listen for and interpret mouse and key events on the configuration display.
     * Holding the "a" key down while performing a mouse button action causes selection of the nearest
     * atom to the cursor and firing of a DisplayPhaseEvent with this atom.
     * Pressing of "s", "b", or "o" keys while display has focus invokes actions that affect the display.
     */
    private class InputEventHandler implements MouseListener, MouseMotionListener, KeyListener,  java.io.Serializable {
        
        Space.Vector point;
        DisplayPhaseEvent dpe;
        
        //not yet configured to do molecule selections
        private boolean atomSelectEnabled = false;
        private boolean moleculeSelectEnabled = false;
        private boolean atomSelected = false;
        private boolean moleculeSelected = false;
        
        InputEventHandler() {
            point = parentSimulation().space().makeVector();
            dpe = new DisplayPhaseEvent(DisplayPhase.this);
        }
        
        public void mouseClicked(MouseEvent evt) {}
        public void mouseEntered(MouseEvent evt) {canvas.requestFocus();}
        public void mouseExited(MouseEvent evt) {canvas.transferFocus();}
        public void mousePressed(MouseEvent evt) {mouseAction(evt);}
        public void mouseReleased(MouseEvent evt) {
            mouseAction(evt);
            dpe.setAtom(null);
            dpe.setMolecule(null);
            atomSelected = false;
            moleculeSelected = false;
        }
        public void mouseDragged(MouseEvent evt) {
            if(atomSelected || moleculeSelected) mouseAction(evt);
            
        }
        public void mouseMoved(MouseEvent evt) {}
        
        private void mouseAction(MouseEvent evt) {
            double toPixels = scale*BaseUnit.Length.Sim.TO_PIXELS;
            double x = (evt.getX() - centralOrigin[0])/toPixels;
            double y = (evt.getY() - centralOrigin[1])/toPixels;
            point.setComponent(0, x);
            point.setComponent(1, y);
            phase().boundary().centralImage(point);
            dpe.setPhase(phase());
            dpe.setPoint(point);
            dpe.setKeyEvent(null);
            dpe.setMouseEvent(evt);
            if(atomSelectEnabled && !atomSelected) {
                dpe.setAtom(selectAtom());
                atomSelected = true;
            }
            if(moleculeSelectEnabled && !moleculeSelected) {
                dpe.setMolecule(selectMolecule());
                moleculeSelected = true;
            }
            fireDisplayPhaseEvent(dpe);
        }
        
        /**
         * Returns the atom nearest the currently selected point
         */
        private Atom selectAtom() {
            Atom nearestAtom = null;
            double r2Min = Double.MAX_VALUE;
            atomIterator.reset();
            while(atomIterator.hasNext()) {
                Atom atom = atomIterator.next();
                double r2 = parentSimulation().space().r2(point,atom.r,phase().boundary());
                if(r2 < r2Min) {
                    nearestAtom = atom;
                    r2Min = r2;
                }
            }
            return nearestAtom;
        }
        
        /**
         * Returns the molecule nearest the currently selected point
         */
        private Molecule selectMolecule() {
            Molecule nearestMolecule = null;
            double r2Min = Double.MAX_VALUE;
            for(Molecule m=phase().firstMolecule(); m!=null; m=m.nextMolecule()) {
                double r2 = parentSimulation().space().r2(point,m.position(),phase().boundary());
                if(r2 < r2Min) {
                   nearestMolecule = m;
                   r2Min = r2;
                }
            }
            return nearestMolecule;
        }  
        
        
        public void keyPressed(KeyEvent evt) {
            char c = evt.getKeyChar();
            if(Character.isDigit(c)) {}
            else if(Character.isLetter(c)) {
                switch(c) {
                    case 'a':
                        atomSelectEnabled = true;
                        moleculeSelectEnabled = false;
                        break;
                    case 'm':
                        atomSelectEnabled = false;
                        moleculeSelectEnabled = true;
                        break;
                    default:
                        break;
                }
            }
            keyAction(evt);
        }
        public void keyReleased(KeyEvent evt) {
            atomSelectEnabled = false;
            moleculeSelectEnabled = false;
            keyAction(evt);
        }
        public void keyTyped(KeyEvent evt) {
            char c = evt.getKeyChar();
            if(Character.isDigit(c)) {setImageShells(Character.getNumericValue(c));}
            else if(Character.isLetter(c)) {
                switch(c) {
                    case 's':
                        canvas.writeScale = !canvas.writeScale;
                        break;
                    case 'o':
                        drawOverflow = !drawOverflow;
                        break;
                    case 'b':
                        setDrawBoundary(!getDrawBoundary());
                        break;
                    default:
                        break;
                }
            }
            keyAction(evt);
        }
        
        private void keyAction(KeyEvent evt) {
            dpe.setPhase(phase());
            dpe.setKeyEvent(evt);
            dpe.setMouseEvent(null);
            fireDisplayPhaseEvent(dpe);
        }
            
    }//end of InputEventHandler
}
