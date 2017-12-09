/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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

import etomica.action.activity.Controller;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.atom.AtomFilter;
import etomica.atom.DiameterHash;
import etomica.atom.DiameterHashByElement;
import etomica.atom.DiameterHashByElementType;
import etomica.space.Space;
import etomica.units.Pixel;

/**
 * Displays a picture of a box, with configurations of molecules, boundaries, and other objects as appropriate, assuming 2-dimensional system.  
 * Instantiates a ConfigurationCanvas (an inner class of this one) for most of the work.
 * DisplayBox is an input event (mouse and key) listener of the canvas.  It receives these 
 * events and uses information from them to form and fire a DisplayBoxEvent to registered listeners.
 *
 * @author David Kofke
 * @author Steve Hotchkiss
 */
 
public class DisplayBox extends Display {
        
    public static final int LEFT = -1;   //Class variables to code for alignment of drawn image within display region
    public static final int CENTER = 0;
    public static final int RIGHT = +1;
    public static final int TOP = -1;
    public static final int BOTTOM = +1;
    //Explicit to 2D because drawing to 2D image
    private final int D = 2;
    protected ColorScheme colorScheme;
    protected DiameterHash diameterHash;
    protected AtomFilter atomFilter = null;
    protected boolean displayBoundary = true;
    LinkedList drawables = new LinkedList();  //was ArrayList before Java2 conversion
    private Box box;
    private boolean graphicResizable = true;
    private final Space space;
    private double sigma = 1.0;
    protected final Simulation sim;
            
    //do not instantiate here; instead must be in graphic method
    public DisplayCanvas canvas = null;

    public final int[] align = new int[D];
    
    /**
     * Size of drawing region of central image, in pixels
     *
     * @see #computeImageParameters()
     */
    protected final int[] drawSize = new int[D];
   
    /**
     * Factor used to scale the size of the image. May be used
     * to scale up or down the image within one box without affecting those
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

    /**
     * When using periodic boundaries, image molecules near the cell boundaries often have parts that overflow
     * into the central cell.  When the box is drawn, these "overflow portions" are not normally
     * included in the central image.  Setting this flag to <code>true</code> causes extra drawing
     * to be done so that the overflow portions are properly rendered.  This is particularly helpful
     * to have on when imageShells is non-zero.  Default value is <code>false</code>.
     */
    private boolean drawOverflow = false;
    
    protected final Controller controller;
  
    /**
     * Warning: after instantiation, clients using G3DSys may need to toggle
     * display.canvas.setVisible false and then true to fix the 'sometimes
     * gray' bug.
     * 
     * i.e.; = new ColorSchemeByType()
     * if(display.canvas instanceof JComponent) {
     * ((JComponent)display.canvas).setVisible(false);
     * ((JComponent)display.canvas).setVisible(true);
     * }
     * @param box
     */
    public DisplayBox(Simulation sim, Box box) {
        super();
        this.sim = sim;
        this.controller = sim.getController();
        this.space = sim.getSpace();
        colorScheme = new ColorSchemeByType(sim);
        setLabel("Configuration");

        align[0] = align[1] = CENTER;

        diameterHash = new DiameterHashByElementType(sim);
        DiameterHashByElement.populateVDWDiameters(((DiameterHashByElementType)diameterHash).getDiameterHashByElement());

        setBox(box);
        setPixelUnit(new Pixel(10));
    }

    /**
     * Sets the size of the box graphic.
     * 
     * @param width  : width of the box graphic
     * @param height : height of the box graphic
     */
    private void setSize(int width, int height) {
    	if(width <= 0 || height <= 0) return;
        java.awt.Dimension temp = new java.awt.Dimension(width, height);
        canvas.setSize(width, height);
        canvas.setMinimumSize(temp);
        canvas.setMaximumSize(temp);
        canvas.setPreferredSize(temp);
        canvas.reshape(width, height);

    }

    /**
     * Gets the size of the box graphic.
     * 
     * @return java.awt.Dimension  : Dimension object containing box graphic size
     */
    public java.awt.Dimension getSize() {
    	return canvas.getSize();
    }

    /**
     * Amount (in pixels) of a simple shift (translation) applied in determing origin.
     * Usually zero, but can be set to any value by passing in dimension and the number
     * of pixels to shift in the specified dimension.
     * <p>
     * Dimension :
     *  0  : X direction
     *  1  : Y direction
     *  2  : Z direction (not used)
     *  
     */
    public void setOriginShift(int dimension, int shift) {
        if(dimension <= D) {
        	originShift[dimension] = shift;
        	java.awt.Dimension dim = canvas.getSize();
        	switch(dimension) {
        	case 0:
        		setSize(dim.width + 2*Math.abs(shift), dim.height);
        		break;
        	case 1:
        		setSize(dim.width, dim.height + 2*Math.abs(shift));
        		break;
        	default:
        		break;
        	}
        }
    }

    /**
     *  
     *  @return int[] : 
     */
    public int[] getOrigin() {
        computeImageParameters();
        return centralOrigin;
    }
 
    /**
     *  
     *  @return int[] : 
     */
    public int[] getDrawSize() {
        computeImageParameters();
        return drawSize;
    }

    /**
     *  
     */
    public void setAlign(int i, int value) {
        align[i] = value;
    }

    /**
     *  
     *  @return int : 
     */
    public int getAlign(int i) {return align[i];}

    /**
     *  
     *  @return final boolean : 
     */
    public final boolean getDrawOverflow() {return drawOverflow;}

    /**
     *  
     */
    public final void setDrawOverflow(boolean b) {drawOverflow = b;}

    /**
     *  
     *  @return double 
     */
    public double getToPixels() {return(canvas.getPixelUnit().toPixels());}

    /**
     *  
     */
    public void setScale(double s) {
        if(s>0) {
            scale = s;
        }
    }

    /**
     *  
     *  @return double 
     */
    public double getScale() {return scale;}

    /**
     *  
     */
    public void addDrawable(Drawable obj) {
        drawables.add(obj);
    }
    /**
     *  
     */
    public void removeDrawable(Drawable obj) {
        drawables.remove(obj);
    }
    /**
     *  
     */
    public void addDrawable(Object obj) {
        if(space.D() == 3) drawables.add(obj);
    }
    /**
     *  
     */
    public void removeDrawable(Object obj) {
        drawables.remove(obj);
    }
    
    /**
     * @return Box : the box associated with this display
     */
    public final Box getBox() {return box;}

    /**
     * Returns the amount of padding added around the edge of the drawing area.
     * This allows atoms at the edge of the boundary to be fully drawn.  The
     * padding is given in terms of simulations length units (Angstroms) and
     * should (in general) be equal to the diameter of atoms in the box.
     */
    public double getPaddingSigma() {
        return sigma;
    }

    /**
     * Sets the amount of padding added around the edge of the drawing area.
     * This allows atoms at the edge of the boundary to be fully drawn.  The
     * padding is given in terms of simulations length units (Angstroms) and
     * should (in general) be equal to the diameter of atoms in the box.
     * DisplayBoxCanvas2D is currently the only thing that pays attention to
     * this.
     */
    public void setPaddingSigma(double sigma) {
        this.sigma = sigma;
    }

    /**
     * Specifies the box for this display.  Updates atomIterator appropriately.
     */
    public void setBox(Box p) {

    	Box oldBox = box;
    	box = p;
    	if(p == null) {
            canvas = null;
            return;
        }

    	double toPixels = (canvas != null ? canvas.getPixelUnit().toPixels() : 10) * scale;
        int boxX = (int)(box.getBoundary().getBoxSize().getX(0) * toPixels + 1);
        int boxY = 1;

        switch(space.D()) {
            case 3:
                boxY = (int)(box.getBoundary().getBoxSize().getX(1) * toPixels);
                if (boxX > boxY) {
                    boxY = boxX;
                }
                else {
                    boxX = boxY;
                }
                boxX *=1.4;
                boxY *=1.4;
                if(canvas == null) {
                	canvas = new DisplayBoxCanvasG3DSys(sim, this, space, controller);
                    setSize(boxX, boxY);
                }
                else {
                	canvas.stopRotate();
                    ((DisplayBoxCanvasG3DSys)canvas).removeObjectByBox(oldBox);
                    ((DisplayBoxCanvasG3DSys)canvas).refreshAtomAgentMgr();
                }
                break;
            case 2:
                boxY = (int)(box.getBoundary().getBoxSize().getX(1) * toPixels + 1);
                canvas = new DisplayBoxCanvas2D(this, space, controller);
                setSize(boxX, boxY);
                break;
            case 1:
            default:
                boxY = drawingHeight;
                canvas = new DisplayBoxCanvas1D(space, this, controller);
                setSize(boxX, boxY);
                break;
        }

        if (graphicResizable == true) {
             setSize(boxX, boxY);
        }

        InputEventHandler listener = new InputEventHandler();
        canvas.addMouseListener(listener);
        canvas.addMouseMotionListener(listener);
        canvas.addKeyListener(listener);

        canvas.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent evt) {
                if((evt.getModifiers() & InputEvent.BUTTON3_MASK) != 0) {
//                    if(DeviceConfigurationEditor.exists) return;
//                    Device editor = new DeviceConfigurationEditor(DisplayBox.this);
//                    ((SimulationGraphic)parentSimulation()).panel().add(editor.graphic(null));
//                    ((SimulationGraphic)parentSimulation()).panel().validate();
//                    ((SimulationGraphic)parentSimulation()).panel().repaint();
                }
            }
        });
    }

    /**
     * 
     */
    public void setBoxCanvas(DisplayCanvas boxCanvas) {
        canvas = boxCanvas;
        if (boxCanvas == null) return;
        if (box == null) throw new IllegalStateException("Cannot set canvas before setting box");
        
        Pixel pixel = canvas.getPixelUnit();
        if (pixel == null) {
            // set a default.  canvases don't set a pixel in their constructors.
            pixel = new Pixel(10);
        }
        int boxX = (int)(box.getBoundary().getBoxSize().getX(0) * pixel.toPixels() * scale + 1);
        int boxY = 1;

        switch(space.D()) {
            case 3:
                boxY = (int)(box.getBoundary().getBoxSize().getX(1) * pixel.toPixels() * scale + 1);
                if (boxX > boxY) {
                    boxY = boxX;
                }
                else {
                    boxX = boxY;
                }
                boxX *=1.4;
                boxY *=1.4;
                break;
            case 2:
                boxY = (int)(box.getBoundary().getBoxSize().getX(1) * pixel.toPixels() * scale + 1);
                break;
            case 1:
            default:
                boxY = drawingHeight;
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
//                    Device editor = new DeviceConfigurationEditor(DisplayBox.this);
//                    ((SimulationGraphic)parentSimulation()).panel().add(editor.graphic(null));
//                    ((SimulationGraphic)parentSimulation()).panel().validate();
//                    ((SimulationGraphic)parentSimulation()).panel().repaint();
                }
            }
        });
    }
    
    /**
     * Mutator method for the color scheme used for this display
     */
    public void setColorScheme(ColorScheme colorScheme) {
        this.colorScheme = colorScheme;
    }
    /**
     * Accessor method for the color scheme used for this display
     * @return ColorScheme
     */
    public ColorScheme getColorScheme() {return colorScheme;}

    /**
     * Returns the DiameterHash used by this DisplayBox.  The default
     * DiameterHash is an instance of DiameterHashByElementType.
     */
    public DiameterHash getDiameterHash() {
        return diameterHash;
    }

    /**
     * Sets the DiameterHash used by this DisplayBox.
     */
    public void setDiameterHash(DiameterHash newDiameterManager) {
        diameterHash = newDiameterManager;
    }

    /**
     * Mutator method for the atom filter that determines which atoms 
     * are displayed.  Atoms for which the filter returns false are not displayed.
     * Default is null, meaning all atoms are displayed.
     */
    public void setAtomFilter(AtomFilter filter) {
        atomFilter = filter;
    }

    /**
     * Accessor method for the atom filter that determines which atoms 
     * are displayed.  Atoms for which the filter returns false are not displayed.
     * Default is null, meaning all atoms are displayed.
     * @return Atomfilter
     */
    public AtomFilter getAtomFilter() {return atomFilter;}

    /**
     *
     * @return LinkedList
     */
    public LinkedList getDrawables() {return(drawables);}
    

    /** 
     * Simulation.GraphicalElement interface method.  Overrides Display method
     * to return the DisplayBox.Canvas as the display object.
     *
     * @param obj ignored by this method.
     * @return Component
     */
    public Component graphic(Object obj) {
        return canvas;
    }

    /**
    * @return int : the current value of imageShells.
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
    
    public void computeImageParameters() {
        int w = canvas.getSize().width - 2 * originShift[0];
        int h = canvas.getSize().height - 2 * originShift[1];
        computeImageParameters2(w, h);
    }

    /**
     * @param w
     * @param h
     */
    public void computeImageParameters2(int w, int h) {
        //Compute factor converting simulation units to pixels for this display
        double toPixels = scale*canvas.getPixelUnit().toPixels();
        //Determine length and width of drawn image, in pixels
        drawSize[0] = (int)(toPixels*getBox().getBoundary().getBoxSize().getX(0));
        drawSize[1] = (space.D()==1) ? drawingHeight: (int)(toPixels*getBox().getBoundary().getBoxSize().getX(1));
        //Find origin for drawing action
        centralOrigin[0] = (int)(getScale()*originShift[0]) + computeOrigin(align[0],drawSize[0],w);
        centralOrigin[1] = (int)(getScale()*originShift[1]) + computeOrigin(align[1],drawSize[1],h);
    }    
      
    private int computeOrigin(int alignX, int drawSizeX, int size) {
        switch(alignX) {
            case   LEFT: return 0;    //same as TOP
            case CENTER: return (size-drawSizeX)/2;
            case  RIGHT: return size-drawSizeX; //same as BOTTOM
            default: return 0;
        }
    }

    /**
     * Repaint the graphic associated with the display box.
     */
    public void repaint() {
        if (canvas != null) {
            canvas.repaint();
        }
    }
      
    //Methods for handling DisplayBoxEvents
    
    /**
     * @return Pixel : unit for conversion between simulation units and display pixels.
     */
    public Pixel getPixelUnit() {
        return canvas.getPixelUnit();
    }

    /**
     * Sets unit for conversion between simulation units and display pixels.
     * @param pixel :
     */
    public void setPixelUnit(Pixel pixel) {
        if(canvas == null) {
            throw new RuntimeException("You gotta at least give me a canvas before you start messing with pixels!");
        }
        canvas.setPixelUnit(pixel);

        int boxX = (int)(box.getBoundary().getBoxSize().getX(0) * pixel.toPixels());
        int boxY = 1;

        switch(space.D()) {
            case 3:
                boxY = (int)(box.getBoundary().getBoxSize().getX(1) * pixel.toPixels());
                if (boxX > boxY) {
                    boxY = boxX;
                }
                else {
                    boxX = boxY;
                }
                boxX *=1.4;
                boxY *=1.4;
                break;
            case 2:
                boxY = (int)(box.getBoundary().getBoxSize().getX(1) * pixel.toPixels());
                break;
            case 1:
            default:
                boxY = drawingHeight;
                break;
        }
        
        setSize(boxX, boxY);

        computeImageParameters();
    }
    
    /**
     * Set the height of the drawing area (only relevant for 1D)
     * @param newDrawingHeight :
     */
    public void setDrawingHeight(int newDrawingHeight) {
        drawingHeight = newDrawingHeight;
    }

    /**
     * return int : height of the drawing area (only relevant for 1D)
     */
    public int getDrawingHeight() {
        return drawingHeight;
    }

    /**
     * Set the flag indicating if the boundary should be drawn.
     * @param b : draw boundary flag
     */
    public void setShowBoundary(boolean b) {
    	displayBoundary = b;
    }

    /**
     * Get the flag indicating if the boundary should be drawn.
     * @return boolean
     */
    public boolean getShowBoundary() {
    	return displayBoundary;
    }

    /**
     * Set the flag indicating whether the graphic should resize from its
     * currently displayed size when a new box is set.
     * @param b
     */
    public void setResizeOnNewBox(boolean b) {
    	graphicResizable = b;
    }

    /**
     * Get the flag indicating whether the graphic should resize from its
     * currently displayed size when a new box is set.
     * @return boolean
     */
    public boolean getResizeOnNewBox() {
    	return graphicResizable;
    }

    /**
     * Number of periodic-image shells to be drawn when drawing this box to the
     * screen.  Default value is 0.
     *
     * @see #repaint
     */
     private int imageShells = 0;
     
     private int drawingHeight = 10;

    /**
     * Class to listen for and interpret mouse and key events on the configuration display.
     * Holding the "a" key down while performing a mouse button action causes selection of the nearest
     * atom to the cursor and firing of a DisplayBoxEvent with this atom.
     * Pressing of "s", "b", or "o" keys while display has focus invokes actions that affect the display.
     */
    private class InputEventHandler implements MouseListener, MouseMotionListener, KeyListener {
        
        Vector point;
        
        //not yet configured to do molecule selections
        private boolean atomSelectEnabled = false;
        private boolean moleculeSelectEnabled = false;
        private boolean atomSelected = false;
        private boolean moleculeSelected = false;
        private boolean rotate = false, zoom = false, translate = false;
        
        InputEventHandler() {
            if(box == null) return;

            point = space.makeVector();
        }
        
        public void mouseClicked(MouseEvent evt) {
            canvas.requestFocus();
        }
        public void mouseEntered(MouseEvent evt) {canvas.requestFocus();}
        public void mouseExited(MouseEvent evt) {canvas.transferFocus();}
        public void mousePressed(MouseEvent evt) {
//			System.out.println("mouse press");
           mouseAction(evt);
            if(space.D() == 3) {
                canvas.setPrevX(evt.getX());
                canvas.setPrevY(evt.getY());
            }
        }
        public void mouseReleased(MouseEvent evt) {
//			System.out.println("mouse release");
            mouseAction(evt);
            atomSelected = false;
            moleculeSelected = false;
        }
         public void mouseDragged(MouseEvent evt) {
//			System.out.println("mouse drag");
            if(atomSelected || moleculeSelected) mouseAction(evt);
           float x = evt.getX();
            float y = evt.getY();
            
            if (rotate  && space.D() == 3) {
                float xtheta = (y - canvas.getPrevY()) * (360f / canvas.getSize().height);
                float ytheta = (x - canvas.getPrevX()) * (360f / canvas.getSize().width);
                canvas.setXRot(canvas.getXRot()+xtheta);
                canvas.setYRot(canvas.getYRot()+ytheta);
            }

            if (translate && space.D() == 3) {
                float xShift = (x - canvas.getPrevX())/-(canvas.getSize().width/canvas.getZoom());
                float yShift = (canvas.getPrevY() - y)/-(canvas.getSize().height/canvas.getZoom());
                canvas.setShiftX(xShift+canvas.getShiftX());
                canvas.setShiftY(yShift+canvas.getShiftY());
            }                                                   

            if (zoom  && space.D() == 3) {
                float xShift = 1f+(x-canvas.getPrevX())/canvas.getSize().width;
                float yShift = 1f+(canvas.getPrevY()-y)/canvas.getSize().height;
                float shift = (xShift+yShift)/2f;
                shift = shift == 1f ? 0: shift < 1f ? shift: -shift;
                canvas.setZoom(canvas.getZoom()+shift);
            }
            
            canvas.repaint();
            if(space.D() == 3) {
                canvas.setPrevX(evt.getX());
                canvas.setPrevY(evt.getY());
            }
            evt.consume();
        }//end of mouseDragged
        
        public void mouseMoved(MouseEvent evt) {}
        
        private void mouseAction(MouseEvent evt) {
            double toPixels = scale * canvas.getPixelUnit().toPixels();
            double x = (evt.getX() - centralOrigin[0])/toPixels;
            double y = (evt.getY() - centralOrigin[1])/toPixels;
            point.setX(0, x);
            point.setX(1, y);
            if(atomSelectEnabled && !atomSelected) {
                atomSelected = true;
            }
            //fire an event if needed
        }
        
        /**
         * Returns the atom nearest the currently selected point
         */
//        private IAtom selectAtom() {
//            IAtom nearestAtom = null;
//            double r2Min = Double.MAX_VALUE;
//            atomIterator.reset();
//            while(atomIterator.hasNext()) {
//                AtomLeaf atom = (AtomLeaf)atomIterator.nextAtom();
//                double r2 = Space.r2(point,atom.getPosition(),getBox().getBoundary());
//                if(r2 < r2Min) {
//                    nearestAtom = atom;
//                    r2Min = r2;
//                }
//            }
//            return nearestAtom;
//        }
        
        /**
         * Returns the molecule nearest the currently selected point
         */
//        private IAtom selectMolecule() {
            //box.moleculeIterator needs to be defined to implement method
//            throw new RuntimeException("method DisplayBox.selectMolecule not implemented");
            
        /*    Atom nearestMolecule = null;
            double r2Min = Double.MAX_VALUE;
            for(AtomIterator iter=box.moleculeIterator; iter.hasNext(); ) {
                Atom m=iter.next();
                double r2 = parentSimulation().space().r2(point,m.coord.position(),box().boundary());
                if(r2 < r2Min) {
                   nearestMolecule = m;
                   r2Min = r2;
                }
            }
            return nearestMolecule;*/
//        }  
        
        
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
                        int drawBoundary = canvas.getDrawBoundary()+1;
                        if (drawBoundary > DisplayCanvas.DRAW_BOUNDARY_MAX) {
                            drawBoundary = 0;
                        }
                        canvas.setDrawBoundary(drawBoundary);
                        break;
                    default:
                        break;
                }
            }
            keyAction(evt);
        }

        private void keyAction(KeyEvent evt) {
            //fire an event if needed
        }
            
    }//end of InputEventHandler
}//end of DisplayBox
