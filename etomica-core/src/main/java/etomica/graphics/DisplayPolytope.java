/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;
import java.awt.Component;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import etomica.action.IAction;
import etomica.space.Vector;
import etomica.math.geometry.Cuboid;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Parallelepiped;
import etomica.math.geometry.PolygonGeneral;
import etomica.math.geometry.Polytope;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.Pixel;

/**
 * Displays a picture of a polytope.  
 * Instantiates a DisplayCanvas for most of the work.
 * DisplayPolytope is an input event (mouse and key) listener of the canvas.
 *
 * @author David Kofke
 * @author Steve Hotchkiss
 */
 
public class DisplayPolytope extends Display implements IAction {
        
    public static final int LEFT = -1;   //Class variables to code for alignment of drawn image within display region
    public static final int CENTER = 0;
    public static final int RIGHT = +1;
    public static final int TOP = -1;
    public static final int BOTTOM = +1;
    public static boolean _3dEnabled;
    private final int D = 2;
    private int drawingHeight = 10;
    private Pixel pixel = new Pixel();
            
    public DisplayCanvas canvas;  //do not instantiate here; instead must be in graphic method

    //Explicit to 2D because drawing to 2D image
    public final int[] align = new int[D];
    
   /**
    * Size of drawing region of central image, in pixels
    *
    */
    protected final int[] drawSize = new int[D];
   
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
      * Factor used to scale the size of the image. May be used
      * to scale up or down the image within one box without affecting those
      * in other displays.  Default value is 1.0.
      */
    protected double scale = 1.0;
          
    private double toPixels;
        
    protected Polytope polytope;
    private final Vector dimensions, maxCoord, minCoord;
  
    public DisplayPolytope(Polytope polytope) {
        super();
        setLabel("Polytope");

        align[0] = align[1] = CENTER;

        dimensions = polytope.getEmbeddedSpace().makeVector();
        minCoord = polytope.getEmbeddedSpace().makeVector();
        maxCoord = polytope.getEmbeddedSpace().makeVector();
        setPolytope(polytope);

        
 //        ((javax.swing.JPanel)graphic()).setLayout(null);
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
     * returns dimensions of polytope as the (max - min) value in
     * each dimension for all the vertices.
     */
    public Vector dimensions() {
        minCoord.E(Double.MAX_VALUE);
        maxCoord.E(-Double.MAX_VALUE);
        final Vector[] vertices = polytope.getVertices();
        for (int i=0; i<vertices.length; i++) {
            for (int j=0; j<dimensions.getD(); j++) {
                if (maxCoord.getX(j) < vertices[i].getX(j)) {
                    maxCoord.setX(j,vertices[i].getX(j));
                }
                else if (minCoord.getX(j) > vertices[i].getX(j)) {
                    minCoord.setX(j,vertices[i].getX(j));
                }
            }
        }
        dimensions.Ev1Mv2(maxCoord,minCoord);
        return dimensions;
    }
    
    public int[] getOrigin() {
        computeImageParameters();
        return centralOrigin;
    }
    
    public int[] getDrawSize() {
        computeImageParameters();
        return drawSize;
    }

    public void setBounds(int x, int y, int width, int height) {
        canvas.setBounds(x,y,width,height);
    }

    public void setAlign(int i, int value) {
        align[i] = value;
    }
    public int getAlign(int i) {return align[i];}

    public double getToPixels() {return(toPixels);}

    public double getScale() {return scale;}
    public void setScale(double s) {
        if(s>0) {
            scale = s;
        }
    }
    
    /**
     * @return the polytope associated with this display
     */
    public final Polytope getPolytope() {return polytope;}

    /**
     * Specifies the polytope for this display.
     */
    public void setPolytope(Polytope p) {
        if(p == null) return;
        polytope = p;

        switch(polytope.getEmbeddedSpace().D()) {
            case 3:
                canvas = new DisplayPolytopeCanvasG3DSys(this);
                break;
            case 2:
                canvas = new DisplayPolytopeCanvas2D(this);
                break;
            case 1:
//                canvas = new DisplayBoundaryCanvas1D(this);
//                break;
            default:
                throw new IllegalArgumentException("can't handle "+polytope.getEmbeddedSpace().D()+" dimensions");
        }
        
        setSize(400, 400);

        InputEventHandler listener = new InputEventHandler();
        canvas.addMouseListener(listener);
        canvas.addMouseMotionListener(listener);
        canvas.addKeyListener(listener);
        
    }

    public void setPolytopeCanvas(DisplayCanvas polytopeCanvas) {
        canvas = polytopeCanvas;
        if (polytopeCanvas == null) return;
        if(polytope == null) throw new IllegalStateException("Cannot set canvas before setting polytope");

        setSize(400, 400);

        InputEventHandler listener = new InputEventHandler();
        canvas.addMouseListener(listener);
        canvas.addMouseMotionListener(listener);
        canvas.addKeyListener(listener);
    }
    
    /** 
     * Simulation.GraphicalElement interface method.  Overrides Display method
     * to return the DisplayBox.Canvas as the display object.
     *
     * @param obj ignored by this method.
     */
    public Component graphic(Object obj) {
        return canvas;
    }

    protected void computeImageParameters() {
        int w = canvas.getSize().width;
        int h = canvas.getSize().height;
        computeImageParameters2(w, h);
    }
    public void computeImageParameters2(int w, int h) {
        //Compute factor converting simulation units to pixels for this display
        toPixels = scale*pixel.toPixels();
        //Determine length and width of drawn image, in pixels
        drawSize[0] = (int)(toPixels*dimensions().getX(0));
        drawSize[1] = (polytope.getEmbeddedSpace().D()==1) ? drawingHeight: (int)(toPixels*dimensions().getX(1));
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

    public void doUpdate() {}
    public void repaint() {canvas.repaint();}
      
    public void setMovable(boolean b) {canvas.setMovable(b);}
    public boolean isMovable() {return canvas.isMovable();}
    public void setResizable(boolean b) {canvas.setResizable(b);}
    public boolean isResizable() {return canvas.isResizable();}
    
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
        if(canvas != null) {
            canvas.setPixelUnit(pixel);
            computeImageParameters();
        }
    }

    /**
     * Class to listen for and interpret mouse and key events on the configuration display.
     */
    private class InputEventHandler implements MouseListener, MouseMotionListener, KeyListener,  java.io.Serializable {
        
        private boolean rotate = false, zoom = false, translate = false;
        
        InputEventHandler() {
            if(polytope == null) return;
        }
        
        public void mouseClicked(MouseEvent evt) {
            canvas.requestFocus();
            //if(parentSimulation().space().D() == 3 && Default.DISPLAY_USE_OPENGL)
            //((DisplayBoxCanvas3DOpenGL)canvas).start();
        }
        public void mouseEntered(MouseEvent evt) {canvas.requestFocus();}
        public void mouseExited(MouseEvent evt) {canvas.transferFocus();}
        public void mousePressed(MouseEvent evt) {
//			System.out.println("mouse press");
            if(polytope.getEmbeddedSpace().D() == 3) {
                canvas.setPrevX(evt.getX());
                canvas.setPrevY(evt.getY());
            }
        }
        public void mouseReleased(MouseEvent evt) {}
         public void mouseDragged(MouseEvent evt) {
//			System.out.println("mouse drag");
           float x = evt.getX();
            float y = evt.getY();
            
            if (rotate  && polytope.getEmbeddedSpace().D() == 3) {
                float xtheta = (y - canvas.getPrevY()) * (360f / canvas.getSize().height);
                float ytheta = (x - canvas.getPrevX()) * (360f / canvas.getSize().width);
                canvas.setXRot(canvas.getXRot()+xtheta);
                canvas.setYRot(canvas.getYRot()+ytheta);
            }

            if (translate && polytope.getEmbeddedSpace().D() == 3) {
                float xShift = (x - canvas.getPrevX())/-(canvas.getSize().width/canvas.getZoom());
                float yShift = (canvas.getPrevY() - y)/-(canvas.getSize().height/canvas.getZoom());
                canvas.setShiftX(xShift+canvas.getShiftX());
                canvas.setShiftY(yShift+canvas.getShiftY());
            }                                                   

            if (zoom  && polytope.getEmbeddedSpace().D() == 3) {
                float xShift = 1f+(x-canvas.getPrevX())/canvas.getSize().width;
                float yShift = 1f+(canvas.getPrevY()-y)/canvas.getSize().height;
                float shift = (xShift+yShift)/2f;
                shift = shift == 1f ? 0: shift < 1f ? shift: -shift;
                canvas.setZoom(canvas.getZoom()+shift);
            }
            
            canvas.repaint();
            if(polytope.getEmbeddedSpace().D() == 3) {
                canvas.setPrevX(evt.getX());
                canvas.setPrevY(evt.getY());
            }
            evt.consume();
        }//end of mouseDragged
        
        public void mouseMoved(MouseEvent evt) {}
        
		public void keyPressed(KeyEvent evt) {
//			System.out.println("key pressed");
			char c = evt.getKeyChar();
			if(Character.isDigit(c)) {}
			else if(Character.isLetter(c)) {
				switch(c) {
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
		}
        public void keyReleased(KeyEvent evt) {}
        public void keyTyped(KeyEvent evt) {}
            
    }//end of InputEventHandler
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
            Space3D space3D = Space3D.getInstance();
            Cuboid cuboid = new Cuboid(space3D);
            Vector3D a = (Vector3D)space3D.makeVector(new double[]{1.0,0.2,0.0});
            Vector3D b = (Vector3D)space3D.makeVector(new double[]{0.0,1.0,0.2});
            Vector3D c = (Vector3D)space3D.makeVector(new double[]{0.2,0.0,1.0});
            Parallelepiped parallelepiped = new Parallelepiped(space3D,a,b,c);
            Space2D space2D = Space2D.getInstance();
            LineSegment lineS1 = new LineSegment(space2D);
            lineS1.setVertex2(space2D.makeVector(new double[]{0,5}));
            LineSegment lineS2 = new LineSegment(space2D);
            lineS2.setVertex1(space2D.makeVector(new double[]{0,5}));
            lineS2.setVertex2(space2D.makeVector(new double[]{10,10}));
            LineSegment lineS3 = new LineSegment(space2D);
            lineS3.setVertex1(space2D.makeVector(new double[]{10,10}));
            lineS3.setVertex2(space2D.makeVector(new double[]{10,0}));
            LineSegment lineS4 = new LineSegment(space2D);
            lineS4.setVertex1(space2D.makeVector(new double[]{0,10}));
            PolygonGeneral quad = new PolygonGeneral(new LineSegment[]{lineS1,lineS2,lineS3,lineS4});
            DisplayPolytope displayPolytope = new DisplayPolytope(cuboid);
            getContentPane().add(displayPolytope.graphic());
        }
    }//end of Applet

    
}//end of DisplayBox
