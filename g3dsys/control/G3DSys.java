package g3dsys.control;

import g3dsys.images.Figure;

import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

import org.jmol.g3d.Graphics3D;

/*
 * index iterator interface in g3dsys
 * etomica classes implement that
 * instead of only vectors, get vectors plus iterator
 * when number of shells change, call setsize on iterator
 */

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
  private FigureManager fm;
  private Display dm;
  private TransformManager tm;

  private java.awt.Container parent; //awt object in which system is contained
  private Graphics3D g3d; //instance of G3D

  private Point3i tempp; //for storing transformations
  
  public G3DSys(java.awt.Container window) {
    //init infrastructure
    parent = window;
    dm = new Display(this);
    dm.setSize(parent.getSize());
    
    //init g3d
    //3/5/2007 rearranged init sequence so g3d is set earlier
    //  needed for imageshell
    g3d = new Graphics3D(dm);
    g3d.setWindowSize(window.getWidth()-10, window.getHeight()-60, false);
    g3d.setBackgroundArgb(0xFF000000);
    g3d.setSlabAndDepthValues(0, Integer.MAX_VALUE);
    
    fm = new FigureManager(this);
    tm = new TransformManager(this);
    tm.clear();
    tm.zoomToPercent(100f);
    tm.setDefaultRotation();

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
    tm.setScreenDimension(parent.getWidth(), parent.getHeight());
    g3d.setWindowSize(parent.getWidth(), parent.getHeight(), false);
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
  /** Gets the pixel depth of the display (zbuffer)
   *  @return depth of the g3d display area in pixels */
  public int getPixelDepth() {
    //used for slab/depth changes; need a better way to do this
    //zoom level affects this
    calcZbufferLimits();
    return zbufferBack - zbufferFront;
  }
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
  private boolean flipYAxis = false;
  private boolean flipXAxis = false;
  public boolean isYFlipped() { return flipYAxis; }
  public boolean isXFlipped() { return flipXAxis; }
  public void setYFlipped(boolean y) { flipYAxis = y; }
  public void setXFlipped(boolean x) { flipXAxis = x; }

  /**
   * Rotate the display around the x and y axes
   * @param degx the number of degrees to rotate x
   * @param degy the number of degrees to rotate y
   */
  public void rotateByXY(float degx, float degy) {
    tm.rotateXYBy(
        (flipXAxis ? -1 : 1)*degx,
        (flipYAxis ? -1 : 1)*degy);
  }
  /**Rotate the display around the x axis
   * @param degrees the number of degrees to rotate
   */
  public void rotateByX(float degrees) {
    rotateByXY(degrees,0);
  }

  /**Rotate the display around the y axis
   * @param degrees the number of degrees to rotate
   */
  public void rotateByY(float degrees) {
    rotateByXY(0,degrees);
  }

  /**Rotate the display around the z axis
   * @param degrees the number of degrees to rotate
   */
  public void rotateByZ(float degrees) { tm.rotateZBy((int)degrees); }

  /** Remove all rotation (and translation) and return to the home position */
  public void rotateToHome() { tm.homePosition(); }
  /** Translate model 10% of pixel width to the left (-x) */
  public void xlateLeft() { tm.translateXYBy(-(int)(tm.width*.1), 0); }
  /** Translate model 10% of pixel width to the right (+x) */
  public void xlateRight() { tm.translateXYBy((int)(tm.height*.1), 0); }
  /** Translate model 10% of pixel width up (-y) */
  public void xlateUp() { tm.translateXYBy(0, -(int)(tm.height*.1)); }
  /** Translate model 10% of pixel width down (+y) */
  public void xlateDown() { tm.translateXYBy(0, (int)(tm.height*.1)); }
  /** Trnaslate model x,y pixels */
  public void xlateXY(int x, int y) { tm.translateXYBy(x, y); }
  /** Set the molecule space location around which rotation occurs */
  public void setCenterOfRotation(Point3f p) {
    //looks like it should work, but haven't tried it yet 2/9/07
    tm.setRotationPointXY(p);
  }
  public Point3f getCenterOfRotation() {
    return tm.getRotationCenter();
  }
  
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
    if(fm.getMinX() == fm.getMaxX()) {
      //initialize; one dimension is flat
      fm.setBoundingBox(minx, miny, minz, maxx, maxy, maxz);
      tm.homePosition();
    }
    else {
      fm.setBoundingBox(minx, miny, minz, maxx, maxy, maxz);
       /* allows the rotation radius to grow along with the bounding box
        * fixes the clipping bug for a growing system, sort of.
        * 
        */
      tm.setDefaultRotation();
    }
  }

  public Point3f getBoundingBoxCenter() {
    return fm.getBoundingBoxCenter();
  }
  
  /* ****************************************************************
   * TransformManager delegation
   * ****************************************************************/
  /**
   * Converts molspace point p (in Angstroms) to a point on the display
   * in pixels.
   * @param p the point to convert
   * @return the pixel coordinate
   */
  public void screenSpace(Point3f p, Point3i s) {
    tempp = tm.transformPoint(p);
    s.set(tempp.x, tempp.y, tempp.z);
  }
  /**
   * Handles perspective based on depth and diameter.
   * @param z depth of an object
   * @param d diameter of an object
   * @return the new perspective-scaled diameter for the object
   */
  public short perspective(int z, float d) {
    return tm.scaleToScreen(z, (int)(d*1000f));
  }
  /**
   * Recalculates the pixel per Angstrom ratio
   */
  public void recalcPPA() {
    //no longer needed since move to TransformManager
  }
  /**
   * Increase the zoom level by i percent
   * @param i the percentage amount to increase 
   */
  public void zoomUp(int i) {
    tm.zoomToPercent(tm.getZoomPercent()+i);
  }
  /**
   * Decrease the zoom level by i percent
   * @param i the percentage amount to decrease
   */
  public void zoomDown(int i) {
    tm.zoomToPercent(tm.getZoomPercentFloat()-i);
  }
  public void setPerspectiveDepth(boolean b) { tm.setPerspectiveDepth(b); }
  public boolean getPerspectiveDepth() { return tm.getPerspectiveDepth(); }



  /* ****************************************************************
   * FigureManager delegation
   * ****************************************************************/
  public float getMinX() { return fm.getMinX(); }
  public float getMinY() { return fm.getMinY(); }
  public float getMinZ() { return fm.getMinZ(); }
  public float getMaxX() { return fm.getMaxX(); }
  public float getMaxY() { return fm.getMaxY(); }
  public float getMaxZ() { return fm.getMaxZ(); }
  public void draw() {
    tm.finalizeTransformParameters();
    g3d.setSlabAndDepthValues((int)tm.slabValue, (int)tm.depthValue);
    fm.draw();
  }

  /**
   * Sets the boundary vectors to be used for calculating image shells
   * @param values the vectors to use
   */
  public void setBoundaryVectors(double[] values) {
    fm.setBoundaryVectors(values);
  }

  /**
   * Sets the index iterator to use for calculating image shells
   * @param i the iterator to use
   */
  public void setBoundaryVectorsIterator(IndexIterator i) {
    fm.setBoundaryVectorsIterator(i);
  }

  /**
   * Gets the index iterator used for calculating image shells
   * @return returns the index iterator
   */
  public IndexIterator getBoundaryVectorsIterator() {
    return fm.getBoundaryVectorsIterator();
  }
  
  /**
   * Gets the boundary vectors used for calculating image shells
   * @return returns the boundary vectors
   */
  public double[] getBoundaryVectors() {
    return fm.getBoundaryVectors();
  }
  
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

  /**
   * Finds the atom distance furthest from the given point; for rotation
   * @param center the reference point
   * @return the furthest distance found, in Angstroms
   */
  public float calcRotationRadius(Point3f center) {
    return fm.calcRotationRadius(center);
  }

  public Point3f getAverageAtomPoint() {
    return fm.getAverageAtomPoint();
  }

  private Point3f origin = new Point3f(0,0,0);
  private int zbufferCenter = 0, zbufferFront = 0, zbufferBack = 0;
  public void calcZbufferLimits() {
    int pixelWidth = (int)(
        java.lang.Math.abs(fm.getMaxX())+
        java.lang.Math.abs(fm.getMinX())*
        tm.scalePixelsPerAngstrom*2);
    zbufferCenter = tm.transformPoint(origin).z;
    zbufferFront = zbufferCenter - pixelWidth;
    zbufferBack = zbufferCenter + pixelWidth;
    System.out.print("zbufferBack: "+zbufferBack);
    System.out.println(", zbufferFront: "+zbufferFront);
  }

  public double getSlabPercent() { return tm.getSlabPercentSetting(); }
  public double getDepthPercent() { return tm.getDepthPercentSetting();  }
  public void setSlabPercent(float val) { tm.slabToPercent(val); }
  public void setDepthPercent(float val) { tm.depthToPercent(val); }

  
  //image shell code follows
  /**
   * For use only by ImageShell class for iteration
   * @return returns the array of current Figures
   */
  public Figure[] getFigs() {
    return fm.getFigs();
  }

  /**
   * Enables or disables image shell
   * @param b boolean value
   */
  public void setEnableImages(boolean b) { fm.setEnableImages(b); }
  /**
   * Returns whether image shell is enabled
   * @return returns enable boolean
   */
  public boolean isEnableImages() { return fm.isEnableImages(); }
  /**
   * Sets the number of image shell layers
   * @param n the number of layers
   */
  public void setLayers(int n) { fm.setLayers(n); }
  /**
   * Gets the number of image shell layers
   * @return returns the number of layers
   */
  public int getLayers() { return fm.getLayers(); }
  /**
   * Sets the boundary drawing style
   * @param b the boundary drawing style to use
   */
  public void setDrawBoundaryType(int b) { fm.setDrawBoundaryType(b); }
  /**
   * Gets the boundary drawing style
   * @return returns the boundary drawing style
   */
  public int getDrawBoundaryType() { return fm.getDrawBoundaryType(); }
  /**
   * Cycles to the next boundary drawing style
   */
  public void cycleDrawBoundaryType() { fm.cycleDrawBoundaryType(); }
  
  public void toggleWireframe() { fm.toggleWireframe(); }

  public void setDefaultRotation() { tm.setDefaultRotation(); }
  
}