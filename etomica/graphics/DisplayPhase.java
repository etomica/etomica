package etomica.graphics;
import java.awt.Component;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.LinkedList;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.Space;
import etomica.action.Action;
import etomica.atom.AtomFilter;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.space.Vector;
import etomica.units.BaseUnit;

/**
 * Displays a picture of a phase, with configurations of molecules, boundaries, and other objects as appropriate, assuming 2-dimensional system.  
 * Instantiates a ConfigurationCanvas (an inner class of this one) for most of the work.
 * DisplayPhase is an input event (mouse and key) listener of the canvas.  It receives these 
 * events and uses information from them to form and fire a DisplayPhaseEvent to registered listeners.
 *
 * @author David Kofke
 * @author Steve Hotchkiss
 */
 
 /* History of changes
  * 09/21/02 (DAK) (temporary) new addDrawable method to add any object, for handling in Space3DOpenGL
  * 01/04/03 (DAK) changed behavior of r/t/z key press to make them "sticky",
  * so that user need not keep key pressed to enable action.  Subsequent press
  * releases key
  * 08/08/03 (DAK) added listener for '<' and '>' keypresses to affect
  * drawExpansionFactor in DisplayPhaseCanvas3DOpenGL
  */
public class DisplayPhase extends Display implements Action, EtomicaElement {
        
    public static final int LEFT = -1;   //Class variables to code for alignment of drawn image within display region
    public static final int CENTER = 0;
    public static final int RIGHT = +1;
    public static final int TOP = -1;
    public static final int BOTTOM = +1;
    public static boolean _3dEnabled;
    private final int D = 2;
    protected ColorScheme colorScheme = new ColorSchemeByType();
    protected AtomFilter atomFilter = AtomFilter.ACCEPT_ALL;
    LinkedList drawables = new LinkedList();  //was ArrayList before Java2 conversion
    private Phase phase;
            
    public DisplayCanvasInterface canvas;  //do not instantiate here; instead must be in graphic method

//Explicit to 2D because drawing to 2D image
    public final int[] align = new int[D];
    
 /**
  * Size of drawing region of central image, in pixels
  *
  * @see #computeDrawSize
  */
    protected final int[] drawSize = new int[D];
   
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
    
    /**
     * Amount of simple shift of drawing origin.
     */
    private final int[] originShift = new int[D];

    private double toPixels;
        
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
   private AtomIterator atomIterator;
   
   static {
 //       _3dEnabled = true;
 //       try {new DisplayPhaseCanvas3DOpenGL(null, 10, 10);}
 //       catch (NoClassDefFoundError err) {_3dEnabled = false;}
   }
  
    public DisplayPhase(Phase phase) {
        super();
        System.out.println("Serenity now");
        setLabel("Configuration");

        align[0] = align[1] = CENTER;

        setPhase(phase);

        
 //        ((javax.swing.JPanel)graphic()).setLayout(null);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Animated display of molecules in a phase as the simulation proceeds");
        return info;
    }
    
    public void setSize(int width, int height) {
        java.awt.Dimension temp = new java.awt.Dimension(width, height);
        canvas.setMinimumSize(temp);
        canvas.setMaximumSize(temp);
        canvas.setPreferredSize(temp);
        canvas.reshape(width, height);
    }
    
    public void actionPerformed() {
        repaint();
    }
    
    
    /**
     * Amount (in pixels) of a simple shift (translation) applied in determing origin.
     * Usually zero, but can be set to any value by setting elements in returned array.
     */
    public int[] getOriginShift() {return originShift;}
        

    public int[] getOrigin() {
        computeImageParameters();
        return centralOrigin;
    }
    
    public int[] getDrawSize() {
        computeImageParameters();
        return drawSize;
    }

    public void setBounds(int x, int y, int width, int height) {
        graphic().setBounds(x,y,width,height);
        canvas.setBounds(x,y,width,height);
    }

    public void setAlign(int i, int value) {
        align[i] = value;
    }
    public int getAlign(int i) {return align[i];}

    public final boolean getDrawOverflow() {return drawOverflow;}
    public final void setDrawOverflow(boolean b) {drawOverflow = b;}

    public double getToPixels() {return(toPixels);}

    public double getScale() {return scale;}
    public void setScale(double s) {
        if(s>0) {
            scale = s;
        }
    }
    
    public void addDrawable(Drawable obj) {
        drawables.add(obj);
    }
    public void removeDrawable(Drawable obj) {
        drawables.remove(obj);
    }
    public void addDrawable(Object obj) {
        if(phase.space().D() == 3) drawables.add(obj);
    }
    public void removeDrawable(Object obj) {
        drawables.remove(obj);
    }
    
    /**
     * @return the phase associated with this display
     */
    public final Phase getPhase() {return phase;}

    /**
     * Specifies the phase for this display.  Updates atomIterator appropriately.
     */
    public void setPhase(Phase p) {
        if(p == null) return;
        phase = p;
        
        int boxX = (int)(phase.boundary().dimensions().x(0) * BaseUnit.Length.Sim.TO_PIXELS);
        int boxY = 1;

        switch(phase.space().D()) {
            case 3:
                boxY = (int)(phase.boundary().dimensions().x(1) * BaseUnit.Length.Sim.TO_PIXELS);
                boxX *=1.4;
                boxY *=1.4;
                    canvas = new DisplayPhaseCanvas3DOpenGL(this, boxX, boxY);
 /*               if(Default.DISPLAY_USE_OPENGL) canvas = new DisplayPhaseCanvas3DOpenGL(this, box, box);
                else canvas = new DisplayPhaseCanvas3DSoftware(this);
 */               break;
            case 2:
                boxY = (int)(phase.boundary().dimensions().x(1) * BaseUnit.Length.Sim.TO_PIXELS);
                canvas = new DisplayPhaseCanvas2D(this);
/*comment this line for applet*/                DefaultGraphic.DISPLAY_USE_OPENGL = false;
                break;
            case 1:
            default:
                canvas = new DisplayPhaseCanvas1D(this);
/*comment this line for applet*/               DefaultGraphic.DISPLAY_USE_OPENGL = false;
                break;
        }
        
        setSize(boxX, boxY);

        InputEventHandler listener = new InputEventHandler();
        canvas.addMouseListener(listener);
        canvas.addMouseMotionListener(listener);
        canvas.addKeyListener(listener);
        
        canvas.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent evt) {
                if((evt.getModifiers() & InputEvent.BUTTON3_MASK) != 0) {
//                    if(DeviceConfigurationEditor.exists) return;
//                    Device editor = new DeviceConfigurationEditor(DisplayPhase.this);
//                    ((SimulationGraphic)parentSimulation()).panel().add(editor.graphic(null));
//                    ((SimulationGraphic)parentSimulation()).panel().validate();
//                    ((SimulationGraphic)parentSimulation()).panel().repaint();
                }
            }
        });

        atomIterator = new AtomIteratorLeafAtoms(p);
    }

    public void setPhaseCanvas(DisplayCanvas phaseCanvas) {
        canvas = phaseCanvas;
        if (phaseCanvas == null) return;
        if(phase == null) throw new IllegalStateException("Cannot set canvas before setting phase");
        
        int boxX = (int)(phase.boundary().dimensions().x(0) * BaseUnit.Length.Sim.TO_PIXELS);
        int boxY = 1;

        switch(phase.space().D()) {
            case 3:
                boxY = (int)(phase.boundary().dimensions().x(1) * BaseUnit.Length.Sim.TO_PIXELS);
                boxX *=1.4;
                boxY *=1.4;
                break;
            case 2:
                boxY = (int)(phase.boundary().dimensions().x(1) * BaseUnit.Length.Sim.TO_PIXELS);
                break;
            case 1:
            default:
                break;
        }
        
        setSize(boxX, boxY);

        InputEventHandler listener = new InputEventHandler();
        canvas.addMouseListener(listener);
        canvas.addMouseMotionListener(listener);
        canvas.addKeyListener(listener);
        
        canvas.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent evt) {
                if((evt.getModifiers() & InputEvent.BUTTON3_MASK) != 0) {
//                    if(DeviceConfigurationEditor.exists) return;
//                    Device editor = new DeviceConfigurationEditor(DisplayPhase.this);
//                    ((SimulationGraphic)parentSimulation()).panel().add(editor.graphic(null));
//                    ((SimulationGraphic)parentSimulation()).panel().validate();
//                    ((SimulationGraphic)parentSimulation()).panel().repaint();
                }
            }
        });
    }
    
    /**
     * Accessor method for the color scheme used for this display
     */
    public void setColorScheme(ColorScheme colorScheme) {
        this.colorScheme = colorScheme;
    }
    /**
     * Accessor method for the color scheme used for this display
     */
    public ColorScheme getColorScheme() {return colorScheme;}
    
    /**
     * Mutator method for the atom filter that determines which atoms 
     * are displayed.  Atoms for which the filter returns false are not displayed.
     * Default is AtomFilter.ALL, according to which all atoms are displayed.
     */
    public void setAtomFilter(AtomFilter filter) {
        atomFilter = (filter == null) ? AtomFilter.ACCEPT_ALL : filter;
        if(canvas != null) canvas.setAtomFilter(atomFilter);
    }
    /**
     * Accessor method for the atom filter that determines which atoms 
     * are displayed.  Atoms for which the filter returns false are not displayed.
     * Default is AtomFilter.ALL, according to which all atoms are displayed.
     */
    public AtomFilter getAtomFilter() {return atomFilter;}
    
    public LinkedList getDrawables() {return(drawables);}
    

    /** 
     * Simulation.GraphicalElement interface method.  Overrides Display method
     * to return the DisplayPhase.Canvas as the display object.
     *
     * @param obj ignored by this method.
     */
    public Component graphic(Object obj) {
        return (Component)canvas;
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
    
/*    public double getFPS() {
        try {
            return Default.DISPLAY_USE_OPENGL ? ((DisplayCanvasOpenGL)canvas).getFps() : 0.;
        }
        catch(NoClassDefFoundError e) {System.out.println("NoClassDefFoundError in getFPS");}
        return 0.0;
    }
    public boolean getUseFpsSleep() {if(Default.DISPLAY_USE_OPENGL)return(((DisplayCanvasOpenGL)canvas).getUseFpsSleep());return(true);}
    public boolean getUseRepaint() {if(Default.DISPLAY_USE_OPENGL)return(((DisplayCanvasOpenGL)canvas).getUseRepaint());return(true);}
    public void setFPS(double fps) {if(Default.DISPLAY_USE_OPENGL)((DisplayCanvasOpenGL)canvas).setAnimateFps(fps);}
    public void setUseFpsSleep(boolean b) {if(Default.DISPLAY_USE_OPENGL)((DisplayCanvasOpenGL)canvas).setUseFpsSleep(b);}
    public void setUseRepaint(boolean b) {if(Default.DISPLAY_USE_OPENGL)((DisplayCanvasOpenGL)canvas).setUseRepaint(b);}
*/
    protected void computeImageParameters() {
        int w = canvas.getSize().width;
        int h = canvas.getSize().height;
        computeImageParameters2(w, h);
    }
    public void computeImageParameters2(int w, int h) {
        //Compute factor converting simulation units to pixels for this display
        toPixels = scale*BaseUnit.Length.Sim.TO_PIXELS;
        //Determine length and width of drawn image, in pixels
        drawSize[0] = (int)(toPixels*getPhase().boundary().dimensions().x(0));
        drawSize[1] = (phase.space().D()==1) ? drawingHeight: (int)(toPixels*getPhase().boundary().dimensions().x(1));
        //Find origin for drawing action
        centralOrigin[0] = (int)(getScale()*originShift[0]) + computeOrigin(align[0],drawSize[0],w);
        centralOrigin[1] = (int)(getScale()*originShift[1]) + computeOrigin(align[1],drawSize[1],h);
    }    
      
    public int computeOrigin(int alignX, int drawSizeX, int size) {
        switch(alignX) {
            case   LEFT: return 0;    //same as TOP
            case CENTER: return (size-drawSizeX)/2;
            case  RIGHT: return size-drawSizeX; //same as BOTTOM
            default: return 0;
        }
    }

    /**
    * Number of periodic-image shells to be drawn when drawing this phase to the
    * screen.  Default value is 0.
    *
    * @see #paint
    */
    private int imageShells = 0;
    
    private int drawingHeight = 10;
     
    public void doUpdate() {}
    public void repaint() {if(!DefaultGraphic.DISPLAY_USE_OPENGL) canvas.repaint();}
      
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
    
    /**
     * Class to listen for and interpret mouse and key events on the configuration display.
     * Holding the "a" key down while performing a mouse button action causes selection of the nearest
     * atom to the cursor and firing of a DisplayPhaseEvent with this atom.
     * Pressing of "s", "b", or "o" keys while display has focus invokes actions that affect the display.
     */
    private class InputEventHandler implements MouseListener, MouseMotionListener, KeyListener,  java.io.Serializable {
        
        Vector point;
        DisplayPhaseEvent dpe;
        
        //not yet configured to do molecule selections
        private boolean atomSelectEnabled = false;
        private boolean moleculeSelectEnabled = false;
        private boolean atomSelected = false;
        private boolean moleculeSelected = false;
        private boolean rotate = false, zoom = false, translate = false;
        
        InputEventHandler() {
            if(phase == null) return;
            point = phase.space().makeVector();
            dpe = new DisplayPhaseEvent(DisplayPhase.this);
        }
        
        public void mouseClicked(MouseEvent evt) {
            canvas.requestFocus();
            //if(parentSimulation().space().D() == 3 && Default.DISPLAY_USE_OPENGL)
            //((DisplayPhaseCanvas3DOpenGL)canvas).start();
        }
        public void mouseEntered(MouseEvent evt) {canvas.requestFocus();}
        public void mouseExited(MouseEvent evt) {canvas.transferFocus();}
        public void mousePressed(MouseEvent evt) {
//			System.out.println("mouse press");
           mouseAction(evt);
            if(phase.space().D() == 3) {
                canvas.setPrevX(evt.getX());
                canvas.setPrevY(evt.getY());
            }
        }
        public void mouseReleased(MouseEvent evt) {
//			System.out.println("mouse release");
            mouseAction(evt);
            dpe.setAtom(null);
            atomSelected = false;
            moleculeSelected = false;
        }
         public void mouseDragged(MouseEvent evt) {
//			System.out.println("mouse drag");
            if(atomSelected || moleculeSelected) mouseAction(evt);
           float x = evt.getX();
            float y = evt.getY();
            
            if (rotate  && phase.space().D() == 3) {
                float xtheta = (y - canvas.getPrevY()) * (360f / canvas.getSize().height);
                float ytheta = (x - canvas.getPrevX()) * (360f / canvas.getSize().width);
                if (!DefaultGraphic.DISPLAY_USE_OPENGL) {
                    ((DisplayPhaseCanvas3DSoftware)canvas).amat.unit();
                    ((DisplayPhaseCanvas3DSoftware)canvas).amat.xrot(xtheta);
                    ((DisplayPhaseCanvas3DSoftware)canvas).amat.yrot(ytheta);
                    ((DisplayPhaseCanvas3DSoftware)canvas).mat.mult(((DisplayPhaseCanvas3DSoftware)canvas).amat);
                } else {
                    canvas.setXRot(canvas.getXRot()+xtheta);
                    canvas.setYRot(canvas.getYRot()+ytheta);
                }
            }

            if (translate && phase.space().D() == 3) {
                float xShift = (x - canvas.getPrevX())/-(canvas.getSize().width/canvas.getZoom());
                float yShift = (canvas.getPrevY() - y)/-(canvas.getSize().height/canvas.getZoom());
                if (!DefaultGraphic.DISPLAY_USE_OPENGL) {
                    ((DisplayPhaseCanvas3DSoftware)canvas).tmat.unit();
                    ((DisplayPhaseCanvas3DSoftware)canvas).tmat.translate(xShift, yShift, 0);
                    ((DisplayPhaseCanvas3DSoftware)canvas).mat.mult(((DisplayPhaseCanvas3DSoftware)canvas).tmat);
               } else {
                    canvas.setShiftX(xShift+canvas.getShiftX());
                    canvas.setShiftY(yShift+canvas.getShiftY());
                }
            }                                                   

            if (zoom  && phase.space().D() == 3) {
                float xShift = 1f+(x-canvas.getPrevX())/canvas.getSize().width;
                float yShift = 1f+(canvas.getPrevY()-y)/canvas.getSize().height;
                float shift = (xShift+yShift)/2f;
                if (!DefaultGraphic.DISPLAY_USE_OPENGL) {
                    ((DisplayPhaseCanvas3DSoftware)canvas).zmat.unit();
                    ((DisplayPhaseCanvas3DSoftware)canvas).zmat.scale(shift, shift, shift);
                    ((DisplayPhaseCanvas3DSoftware)canvas).mat.mult(((DisplayPhaseCanvas3DSoftware)canvas).zmat);
                    setScale(shift*shift);
                    //xfac *= shift*shift;
                    //settings.setAtomScreenScale(xfac);
                    //settings.setBondScreenScale(xfac);
                    //settings.setVectorScreenScale(xfac);
                } else {
                    shift = shift == 1f ? 0: shift < 1f ? shift: -shift;
                    canvas.setZoom(canvas.getZoom()+shift);
                }
            }
            
            if (!DefaultGraphic.DISPLAY_USE_OPENGL) canvas.repaint();
            if(phase.space().D() == 3) {
                canvas.setPrevX(evt.getX());
                canvas.setPrevY(evt.getY());
            }
            evt.consume();
        }//end of mouseDragged
        
        public void mouseMoved(MouseEvent evt) {}
        
        private void mouseAction(MouseEvent evt) {
            double x = (evt.getX() - centralOrigin[0])/toPixels;
            double y = (evt.getY() - centralOrigin[1])/toPixels;
            point.setX(0, x);
            point.setX(1, y);
   //         phase().boundary().centralImage(point);
            dpe.setPhase(getPhase());
            dpe.setPoint(point);
            dpe.setKeyEvent(null);
            dpe.setMouseEvent(evt);
            if(atomSelectEnabled && !atomSelected) {
                dpe.setAtom(selectAtom());
                atomSelected = true;
            }
/*            if(moleculeSelectEnabled && !moleculeSelected) {
                dpe.setMolecule(selectMolecule());
                moleculeSelected = true;
            }*/
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
                Atom atom = atomIterator.nextAtom();
                double r2 = Space.r2(point,atom.coord.position(),getPhase().boundary());
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
        private Atom selectMolecule() {
            //phase.moleculeIterator needs to be defined to implement method
            throw new RuntimeException("method DisplayPhase.selectMolecule not implemented");
            
        /*    Atom nearestMolecule = null;
            double r2Min = Double.MAX_VALUE;
            for(AtomIterator iter=phase.moleculeIterator; iter.hasNext(); ) {
                Atom m=iter.next();
                double r2 = parentSimulation().space().r2(point,m.coord.position(),phase().boundary());
                if(r2 < r2Min) {
                   nearestMolecule = m;
                   r2Min = r2;
                }
            }
            return nearestMolecule;*/
        }  
        
        
//		public void keyPressed(KeyEvent evt) {
//			System.out.println("key pressed");
//			char c = evt.getKeyChar();
//			if(Character.isDigit(c)) {}
//			else if(Character.isLetter(c)) {
//				switch(c) {
//					case 'a':
//						atomSelectEnabled = true;
//						moleculeSelectEnabled = false;
//						break;
//					case 'm':
//						atomSelectEnabled = false;
//						moleculeSelectEnabled = true;
//						break;
//					case 'r':
//						rotate = true;
//						zoom = false;
//						translate = false;
//						break;
//					case 'z':
//						rotate = false;
//						zoom = true;
//						translate = false;
//						break;
//					case 't':
//						rotate = false;
//						zoom = false;
//						translate = true;
//						break;
//				   default:
//					   break;
//				}//end switch
//			}
//			keyAction(evt);
//		}
		public void keyPressed(KeyEvent evt) {
//			System.out.println("key pressed");
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
					case 'r':
						rotate = !rotate;
						zoom = false;
						translate = false;
						break;
					case 'z':
						rotate = false;
						zoom = !zoom;
						translate = false;
						break;
					case 't':
						rotate = false;
						zoom = false;
						translate = !translate;
						break;
				   default:
					   break;
				}//end switch
			}
			keyAction(evt);
		}
        public void keyReleased(KeyEvent evt) {
//        	System.out.println("released");
            atomSelectEnabled = false;
            moleculeSelectEnabled = false;
//            rotate = false;
//            zoom = false;
//            translate = false;
            keyAction(evt);
        }
        public void keyTyped(KeyEvent evt) {
            char c = evt.getKeyChar();
            if(Character.isDigit(c)) {setImageShells(Character.getNumericValue(c));}
            else if(Character.isLetter(c)) {
                switch(c) {
                    case 's':
                        canvas.setWriteScale(!canvas.getWriteScale());
                        break;
                    case 'o':
                        drawOverflow = !drawOverflow;
                        break;
                    case 'b':
                        canvas.setDrawBoundary(canvas.getDrawBoundary()+1);
                        break;
                    case 'q':
                        canvas.setQuality(canvas.getQuality()+1);
//                        canvas.setHighQuality(!canvas.getHighQuality());
                        break;
                    default:
                        break;
                }
            }
            else if(canvas instanceof DisplayPhaseCanvas3DOpenGL){
            	DisplayPhaseCanvas3DOpenGL canvasGL = (DisplayPhaseCanvas3DOpenGL)canvas;
            	switch(c) {
            		case '<':
//						setImageShells(0);
            			canvasGL.setDrawExpansionFactor(canvasGL.getDrawExpansionFactor()-0.02);
            			drawOverflow = false;
            			break;
					case '>':
//						setImageShells(0);
						canvasGL.setDrawExpansionFactor(canvasGL.getDrawExpansionFactor()+0.02);
						drawOverflow = false;
						break;
            	}            			
            }
            keyAction(evt);
        }
        
        private void keyAction(KeyEvent evt) {
            dpe.setPhase(getPhase());
            dpe.setKeyEvent(evt);
            dpe.setMouseEvent(null);
            fireDisplayPhaseEvent(dpe);
        }
            
    }//end of InputEventHandler
}//end of DisplayPhase
