package g3dsys.control;

import g3dsys.images.*;

import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

import org.jmol.g3d.Graphics3D;

/**
 * Abstraction and delegation container class to facilitate communication
 * between the CoordMapper, FigureManager, and G3DPanel without deeply
 * coupling them in their own definitions.
 * 
 * Users need only instantiate G3DSys with the target AWT component and then
 * start adding Figures for them to appear. No Figure constructors should be
 * called directly, as the addFig method will do that itself based on the
 * figure type argument (e.g.; BALL, LINE, BOX).
 */

public class G3DSys {

  /** implemented shapes */
  public static final int BALL = 0;
  public static final int CYLINDER = 1;
  public static final int RECTANGLE = 2;
  public static final int LINE = 3;
  public static final int BOX = 4;
  public static final int AXES = 5;
  public static final int TIMEDBALL = 6;

  // classes needing communication amongst themselves
  // coupled in function not in code
  private CoordMapper cm;
  private FigureManager fm;
  private Display dm;

  private java.awt.Container parent; //awt object in which system is contained
  private Graphics3D g3d; //instance of G3D

  public G3DSys(java.awt.Container window) {
    //init infrastructure
    parent = window;
    dm = new Display(this);

    fm = new FigureManager(this);
    cm = new CoordMapper(this);
    dm.setSize(parent.getSize());

    //init g3d
    g3d = new Graphics3D(dm);
    g3d.setWindowSize(window.getWidth()-10, window.getHeight()-60, false);
    g3d.setBackgroundArgb(0xFF000000);
    g3d.setSlabAndDepthValues(0, Short.MAX_VALUE);

    //enable scaling resize
    parent.addComponentListener(new ComponentListener() {
      public void componentHidden(ComponentEvent e) {}
      public void componentMoved(ComponentEvent e) { fastRefresh(); }
      public void componentShown(ComponentEvent e) {}
      public void componentResized(ComponentEvent e) { refresh(); }});

    //add keyboard controls
    dm.addKeyListener(new KeyboardControl(this));
    //add mouse controls
    new MouseManager(dm,this);
    //dm.addMouseListener(new MouseControl(this));

    parent.add(dm);
    refresh();
    dm.repaint();
    parent.repaint();
  }




  /* ****************************************************************
   * G3D-related methods
   * ****************************************************************/

  /** For a thorough redraw of the display */
  public void refresh() {
    dm.setSize(parent.getSize());
    dm.setLocation(0,0);
    g3d.setWindowSize(parent.getWidth(), parent.getHeight(), false);
    cm.recalcPPA();
    dm.repaint();
  }

  /** A less thorough redraw for a quick refresh */
  public void fastRefresh() {
    //dm.paint(dm.getGraphics());
    dm.repaint();
  }

  /**Change the background color of the G3D display
   * @param color the G3d background color to use
   */
  public void setBGColor(int color) { g3d.setBackgroundArgb(color); }

  /** Get G3D
   *  @return associated Graphics3D object */
  public Graphics3D getG3D() { return g3d; }




  /* ****************************************************************
   * Dimension accessors for calculation pixel per Angstrom ratio
   * ****************************************************************/
  /** Gets the pixel width of the display
   *  @return width of the g3d display area in pixels */
  public int getPixelWidth() { return parent.getWidth(); }
  /** Gets the pixel height of the display
   *  @return height of the g3d display area in pixels */
  public int getPixelHeight() { return parent.getHeight(); }
  /** Gets the pixel depth of the display
   *  @return depth of the g3d display area in pixels */
  public int getPixelDepth() { return cm.angToPixel(fm.getDepth()); }
  /** Gets the Angstrom width of the model
   *  @return width of the model in Angstroms */
  public float getAngstromWidth() { return fm.getWidth(); }
  /** Gets the Angstrom height of the model
   *  @return height of the model in Angstroms */
  public float getAngstromHeight() { return fm.getHeight(); }
  /** Gets the Angstrom depth of the model
   *  @return depth of the model in Angstroms */
  public float getAngstromDepth() { return fm.getDepth(); }




  /* ****************************************************************
   * Methods for modifying the model
   * ****************************************************************/

  /**
   * Adds the Figure to the associated FigureManager
   * @param f the Figure to add
   */
  public void addFig(Figure f) {
    fm.addFig(f);
  }
  
  /**
   * Removes the Figure to from associated FigureManager
   * @param f the Figure to remove
   */
  public Figure removeFig(Figure f) {
    return fm.removeFig(f);
  }
  
  /* ****************************************************************
   * Rotation and translation delegation
   * ****************************************************************/
  private boolean flipYAxis = true;
  private boolean flipXAxis = false;
  /**Rotate the display around the x axis
   * @param degrees the number of degrees to rotate
   */
  public void rotateByX(float degrees) {
    cm.rotateByX( (flipXAxis ? -1 : 1) * degrees);
  }

  /**Rotate the display around the y axis
   * @param degrees the number of degrees to rotate
   */
  public void rotateByY(float degrees) {
    cm.rotateByY( (flipYAxis ? -1 : 1) * degrees);
  }

  /**Rotate the display around the z axis
   * @param degrees the number of degrees to rotate
   */
  public void rotateByZ(float degrees) { cm.rotateByZ(degrees);	}

  /** Remove all rotation (and translation) and return to the home position */
  public void rotateToHome() { cm.rotateToHome(); }
  /** Translate model 20% of pixel width to the left (-x) */
  public void xlateLeft() { cm.xlateX(-1*getPixelWidth()/5); }
  /** Translate model 20% of pixel width to the right (+x) */
  public void xlateRight() { cm.xlateX(getPixelWidth()/5); }
  /** Translate model 20% of pixel width up (-y) */
  public void xlateUp() { cm.xlateY(-1*getPixelHeight()/5); }
  /** Translate model 20% of pixel width down (+y) */
  public void xlateDown() { cm.xlateY(getPixelHeight()/5); }
  /** Translate model horizontally x pixels */
  public void xlateX(int x) { cm.xlateX(x); }
  /** Translate model vertically y pixels */
  public void xlateY(int y) { cm.xlateY(y); }
  /** Set the molecule space location around which rotation occurs */
  public void setCenterOfRotation(Point3f p) { cm.setCenterOfRotation(p); }
  public Point3f getCenterOfRotation() { return cm.getCenterOfRotation(); }


  /**
   * Set the min and maximum coordinates for the model.  The zoom level is set
   * such that everything within the box is viewable.
   * @param minx minimum x
   * @param miny minimum y
   * @param minz minimum z
   * @param maxx maximum x
   * @param maxy maximum y
   * @param maxz maximum z
   */
  public void setBoundingBox(float minx, float miny, float minz,
          float maxx, float maxy, float maxz) {
      fm.setBoundingBox(minx, miny, minz, maxx, maxy, maxz);
  }

  /* ****************************************************************
   * CoordMapper delegation
   * ****************************************************************/
  /**
   * Converts molspace point p (in Angstroms) to a point on the display
   * in pixels.
   * @param p the point to convert
   * @return the pixel coordinate
   */
  public void screenSpace(Point3f p, Point3i s) { cm.screenSpace(p, s); }
  /**
   * Handles perspective based on depth and diameter.
   * @param z depth of an object
   * @param d diameter of an object
   * @return the new perspective-scaled diameter for the object
   */
  public int perspective(int z, float d) { return cm.perspective(z,d); }
  /**
   * Recalculates the pixel per Angstrom ratio
   */
  public void recalcPPA() { cm.recalcPPA(); }
  /**
   * Increase the zoom level by i percent
   * @param i the percentage amount to increase 
   */
  public void zoomUp(int i) { cm.zoomUp(i); }
  /**
   * Decrease the zoom level by i percent
   * @param i the percentage amount to decrease
   */
  public void zoomDown(int i) { cm.zoomDown(i); }




  /* ****************************************************************
   * FigureManager delegation
   * ****************************************************************/
  public float getMinX() { return fm.getMinX(); }
  public float getMinY() { return fm.getMinY(); }
  public float getMinZ() { return fm.getMinZ(); }
  public float getMaxX() { return fm.getMaxX(); }
  public float getMaxY() { return fm.getMaxY(); }
  public float getMaxZ() { return fm.getMaxZ(); }
  public void draw() { fm.draw();}


  /**
   * Converts the given AWT color to an argb int
   * @param color the color to be converted
   * @returns an argb int
   **/
  private static int color2argb(java.awt.Color color) {
    float[] compArray = color.getComponents(null);
    int a = (int)(compArray[3]*255+0.5);
    int r = (int)(compArray[0]*255+0.5);
    int g = (int)(compArray[1]*255+0.5);
    int b = (int)(compArray[2]*255+0.5);
    int argb = (a << 24) | (r << 16) | (g << 8) | b;
    return argb;
  }

  /**
   * Converts the given AWT color to an argb int
   * @param color the color to be converted
   * @returns an argb int
   **/
  public static short getColix(java.awt.Color color) {
    return Graphics3D.getColix(color2argb(color));
  }
}