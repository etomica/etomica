package simulate;
import java.io.*;
import java.awt.*;  //subclasses Component; for Graphics and Color in draw method
import java.beans.*;
import java.util.*;

/**
 * Each Phase must contain exactly one Space.  The dimensionality of the system 
 * is defined here.  The space object prescribes how
 * distances are measured, and it sets bounds on the allowable values of 
 * particle coordinates, defining what happens to molecules that go beyond
 * these bounds.  Thus subclasses of space can be defined to set up (for example)
 * periodic boundary conditions.
 *
 * Space also defines a large set of static utility methods to perform vector
 * arithmetic.  They are written specifically for vectors in a 2-dimensional
 * space.  Their purpose is to permit addition, multiplication, dot-products,
 * etc. of scalars and vectors without using a for-loop construct.  Extension
 * to simulation of 3-dimensional (or other) systems can be accomplished by
 * defining another Space routine with these methods appropriately redefined.
 * Since these are invoked as static methods, this object must replace (rather
 * than extend) the Space object defined here.  Methods used to compute vector
 * differences are written as instance (non-static) methods so they can
 * simply be overridden to re-define how distances are computed.  It is 
 * thought by us that the awkwardness of this approach is offset by the
 * speedup gained in using static, unrolled-loop methods to do these often-used
 * calculations.
 *
 * The utility methods are named according to a code that describes their function.
 *
 *      In the argument list, "u" represents the "left-hand side" if the 
 *      method calculation results in a vector; the array passed in this
 *      position will have its values updated according to the calculation.
 *      If the calculation results in a scalar, this is the return value
 *      of the method and no "u" is present in the argument list.
 *      
 *      v1, v2, etc. represent vectors appearing on the right-hand side
 *      of the vector equation, while a1, a2, etc. represent scalars 
 *      appearing there.
 *      
 *      Method names are coded to suggest how the vector equation
 *      would be written using standard notation.  Thus:
 *      
 *      uE --> "u equals"; returns a vector in the first argument
 *       P --> "plus"
 *       M --> "minus"
 *       T --> "times"
 *       D --> "dot"; inner product
 *       S --> "squared"; vector square magnitude
 *      PE --> "plus equals"  (+=); increments u
 *      ME --> "minus equals" (-=)
 *      underscores represent parentheses ("left paren" omitted if would be first character of method name)
 *
 * Space and its subclasses are Beans, and are directly accessed at design time.
 * At design time, a Space draws itself as a small yellow rectangle.  At run time it
 * may be invisible or it may be set to draw its boundaries.
 *
 * @see Phase
 */

public class Space extends Component {

 /**
  * Dimensionality of the space.  If this value is changed, all vector utility methods
  * must be updated accordingly.
  */
    public static final int D = 2;
 
 /**
  * True for spaces that have periodic boundaries
  */
    public boolean periodic;
  
 
 /**
  * Size of Phase (width, height) in Angstroms
  * Default value is 1.0 for each dimension.
  */
    protected final double[] dimensions = new double[D];
    
 /**
  * Size of Space, in pixels
  *
  * @see #computeDrawSize
  */
    protected final int[] drawSize = new int[D];
  
 /**
  * Factor used to scale the size of the image within the phase. May be used
  * to scale up or down the image within one phase without affecting those
  * in other phases.  Default value is 1.0.
  *
  * @see Phase#paint
  */
    protected double scale = 1.0;
    
    private static final double[] work = new double[D];       //temporary work vector
 
 /**
  * @see SpacePeriodicCubic#getOverflowShifts
  */
    protected final double[][] shift0 = new double[0][D];
   
   /**
    * Phase in which this Space resides
    */
    protected Phase parentPhase;
   
   /**
    * Volume of the phase, in Angstroms^D.
    *
    * @see #computeVolume
    */
    public double volume;
   
   /**
    * Number of period-image shells to be drawn about central image
    */
    protected int nShells = 0;
   
   /**
    * Number of period images to be drawn about the central image.  Determined
    * by the shape/dimensionality of the central image (i.e., the coordination
    * number of the periodic lattice) and the value of nShells.
    */
    protected int nImages;
   
   /**
    * Coordinate origin for central image
    */
    protected final int[] centralOrigin = new int[D];
    
   /**
    * Coordinate origins for all images about central image.
    */
    protected int[][] origins;
    
   /**
    * Creates Space with unit dimensions and computes its volume.
    */
    public Space() {
        dimensions[0] = 1.0;
        dimensions[1] = 1.0;
        computeVolume();
        periodic = false;
    }
    
   /**
    * Sets parentPhase and invokes setDimensions to permit Phase to adjust
    * TO_PIXELS.
    *
    * @param p the phase to be identified as this Space's parentPhase
    * @see #setDimensions
    */
    public void setParentPhase(Phase p) {
        parentPhase = p;
        setDimensions(0,dimensions[0]);
    }
    
    public final boolean isPeriodic() {return periodic;}  
    
   /**
    * Draws a light gray outline of the space if <code>visible</code> is set to
    * <code>true</code>.  The size of the outline is determined by 
    * <code>drawSize[]</code>.
    *
    * @param g      the graphics object to which the image is drawn
    * @param origin the coordinate origin (in pixels) for drawing the image
    * @see Phase#paint
    * @see #computeDrawSize
    */
    public void draw(Graphics g, int[] origin) {
        if(isVisible()) {
            g.setColor(Color.gray.brighter());
            g.drawRect(origin[0],origin[1],drawSize[0]-1,drawSize[1]-1);
        }
    }
    
   /**
    * @return the current value of scale
    * @see #setScale
    */
    public double getScale() {return scale;}
    
   /**
    * Sets scale and nShells and updates drawSize.  Likely to be overridden
    * in subclasses of Space.
    *
    * @param s   the nominal drawing scale.  No action is taken if s <= 0.0.
    * @param n   the new value of nShells.  No action is taken if n < 0.
    * @see #computeDrawSize
    * @see Phase#paint
    */
    public void setScale(double s, int n) {   //likely to override
        if(s>0 && n>=0) {
          scale = s;
          nShells = n;
          computeDrawSize();
        }
    }
   
   /**
    * Updates drawSize according to the current values of scale, dimensions,
    * and Phase.TO_PIXELS.
    * drawsize[i] = Phase.TO_PIXELS * scale * dimensions[i]
    */
    protected void computeDrawSize() {
        drawSize[0] = Phase.toPixels(scale*dimensions[0]);
        drawSize[1] = Phase.toPixels(scale*dimensions[1]);
        resetOrigins(nShells);
    }
    
   /**
    * @return the drawSize[] array
    */
    public int[] getDrawSize() {return drawSize;}

   /**
    * @return the dimensions[] array
    */
    public double[] getDimensions() {return dimensions;}
   
   /**
    * @return the value of dimensions[i]
    * @param i  index of desired value of dimensions[]
    */
    public double getDimensions(int i) {return dimensions[i];}
    
   /**
    * Sets one value in the dimensions array and updates drawSize array
    * and <code>volume</code>
    *
    * @param i   index of the value to be set
    * @param dim the new value of dimensions[i]
    */
    public void setDimensions(int i, double dim) {
        dimensions[i] = dim;
        if(parentPhase != null) {parentPhase.resetTO_PIXELS();}
        computeDrawSize();
        computeVolume();
    }
   
   /**
    * Scales all dimensions by a constant multiplicative factor
    * Calls setDimensions(i,scale*dimensions[i]) for all i=1..D
    *
    * @param scale the scaling factor. 
    */
    public void inflate(double scale) {
        setDimensions(0,scale*dimensions[0]);
        setDimensions(1,scale*dimensions[1]);
    }
   
   /**
    * Computes the volume of the space based on the current values in
    * dimensions[], and assigns the result to <code>volume</code>.
    * Likely to override in subclasses.
    */
    protected void computeVolume() {
        volume = dimensions[0]*dimensions[1];
    }
    
   /**
    * Computes origins needed to draw central and periodic images.
    */
    protected void resetOrigins(int n) {    //likely to override
        nShells = n;
        if(parentPhase != null) resetCentralOrigin(parentPhase.getPhaseSize());
        resetImageOrigins();
    }
    
   /**Likely to override
    * in subclasses.
    */
    public void resetImageOrigins() {
        nImages = 0;
        origins = new int[nImages][];
    }
   
    protected void resetCentralOrigin(int[] phaseSize) {
        Space.uEa1T_v1Mv2_(centralOrigin,0.5,phaseSize,drawSize);  // Maybe get this from space;  origin = 0.5*(phaseSize - spaceSize)
    }
        
    public final int[] getCentralOrigin() {
        return centralOrigin;
    }
    
    public int[] getCopyOrigin() {  //Origin for original copy used to produce periodic images
        return centralOrigin;
    }
    
   /**
    * Takes value for nShells, checks to see if it matches current value,
    * resets value and origins if they differ, then returns imageOrigins array for
    * new nShells value
    *
    * @param n  new value of nShells
    * @return origins[][] array
    * @see Phase#paint
    */
    public int[][] getImageOrigins(int n) {
        if(n != nShells) {resetOrigins(n);}
        return origins;
    }
    
   /**
    * @return shift0
    * @see #shift0
    * @see Phase#paint
    */
    public double[][] getOverflowShifts(double[] r, double distance) {return shift0;}  //called only if periodic
    
   /**
    * Design-time manipulation of dimensions[0].  Used in lieu of an array editor.
    */
    public void setXDimension(double dim) {setDimensions(0, dim);} //used for property list editor
    public double getXDimension() {return dimensions[0];}
    public void setYDimension(double dim) {setDimensions(1,dim);}
    public double getYDimension() {return dimensions[0];}
   
   /**
    * Method to handle placement of molecules that go beyond boundaries of Space.
    * Used by periodic-boundary subclasses.
    * 
    * @see Phase#paint
    */
    public void repositionMolecules(Species species) {return;}  

   /**
    * Method to draw Space at design time.  Draws a yellow rectangle to
    * confirm placement of space in phase.
    *
    * @param g graphic object to which image is drawn.
    */
    public void paint(Graphics g) {
        if(Beans.isDesignTime()) {
            g.setColor(Color.yellow);
            g.fillRect(0,0,getSize().width,getSize().height);
        }
    }
 
    /* The remainder of the methods do vector arithmetic.
    
    */  
    
   /**
    * Assigns the constant a1 to all elements of the array u (u = a1).
    */
    public static void uEa1(double[] u, double a1) {
        u[0] = u[1] = a1;
    }
    
    public static void uEa1(int[] u, int a1) {
        u[0] = u[1] = a1;
    }
    
   /**
    * Element-by-element copy of array v1 to array u (u = v1)
    */
    public static void uEv1(double[] u, double[] v1) {
        u[0] = v1[0];
        u[1] = v1[1];
    }

   /**
    * Element-by-element copy of array v1 to array u (u = v1)
    */
    public static void uEv1(int[] u, int[] v1) {
        u[0] = v1[0];
        u[1] = v1[1];
    }
    
   /**
    * @return the scalar dot product of arrays v1 and v2
    */
    public static double v1Dv2(double[] v1, double[] v2) {
        return (v1[0]*v2[0] + v1[1]*v2[1]);
    }
   
   /**
    * @return the vector squared, the magnitude of the vector v1: v1^2 = v1 . v1
    */
    public static double v1S(double[] v1) {
        return (v1[0]*v1[0] + v1[1]*v1[1]);
    }
   
   /**
    * @return the scalar |v1-v2|^2
    */
    public static double v1Mv2_S(double[] v1, double[] v2) {
        double w0 = v1[0] - v2[0];
        double w1 = v1[1] - v2[1];
        return w0*w0 + w1*w1;
    }        
    
    public static void uEa1Tv1(double[] u, double a1, double[] v1) {
        u[0] = a1*v1[0];
        u[1] = a1*v1[1];
        return;
    }
    
    public static void uEv1Mv2(double[] u, double[] v1, double[] v2) {
        u[0] = v1[0] - v2[0];
        u[1] = v1[1] - v2[1];
        return;
    }
    
    public void uEr1Mr2(double[] u, double[] v1, double[] v2) {  //for overriding by subclasses
        uEv1Mv2(u, v1, v2);
    }
    public double r1Mr2_S(double[] v1, double[] v2) {  //for overriding by subclasses
        return v1Mv2_S(v1, v2);
    }
    public double r1iMr2i(int i, double[] v1, double[] v2) {  //for overriding
        return v1[i] - v2[i];
    }
    
    public static void uPEv1(double[] u, double[] v1) {
        u[0] += v1[0];
        u[1] += v1[1];
        return;
    }
        
    public static void uMEv1(double[] u, double[] v1) {
        u[0] -= v1[0];
        u[1] -= v1[1];
        return;
    }
    
    public static void uPEa1(double[] u, double a1) {
        u[0] += a1;
        u[1] += a1;
    }
        
    public static void uMEa1(double[] u, double a1) {
        u[0] -= a1;
        u[1] -= a1;
    }
        
    public static void uTEa1(double[] u, double a1) {
        u[0] *= a1;
        u[1] *= a1;
        return;
    }
    
    public static void uDEa1(double[] u, double a1) {
        u[0] /= a1;
        u[1] /= a1;
    }
        
    public static void uPEa1Tv1(double[] u, double a1, double[] v1) {
        u[0] += a1*v1[0];
        u[1] += a1*v1[1];
        return;
    }
    
    public static void uMEa1Tv1(double[] u, double a1, double[] v1) {
        u[0] -= a1*v1[0];
        u[1] -= a1*v1[1];
        return;
    }
 
    public static void uEa1Tv1Ma2Tv2(double[] u, double a1, double[] v1, double a2, double[] v2) {
        u[0] = a1*v1[0] - a2*v2[0];
        u[1] = a1*v1[1] - a2*v2[1];
        return;
    }
    
    public static void uEa1Tv1Pa2Tv2(double[] u, double a1, double[] v1, double a2, double[] v2) {
        u[0] = a1*v1[0] + a2*v2[0];
        u[1] = a1*v1[1] + a2*v2[1];
        return;
    }
    
    public static void uEa1T_v1Mv2_(double[] u, double a1, double[] v1, double[] v2) {
        u[0] = a1*(v1[0]-v2[0]);
        u[1] = a1*(v1[1]-v2[1]);
    }
    
    public static void uEa1T_v1Mv2_(int[] u, double a1, int[] v1, int[] v2) {
        u[0] = (int)(a1*(v1[0]-v2[0]));
        u[1] = (int)(a1*(v1[1]-v2[1]));
    }
    
    public static void uEv1Pa1Tv2(double[] u, double[] v1, double a1, double[] v2) {
        u[0] = v1[0] + a1*v2[0];
        u[1] = v1[1] + a1*v2[1];
    }
    
    public static void uEv1Ma1Tv2(double[] u, double[] v1, double a1, double[] v2) {
        u[0] = v1[0] - a1*v2[0];
        u[1] = v1[1] - a1*v2[1];
    }
    
    public static void unitVector(double[] u, double[] v1, double[] v2) {
        uEv1Mv2(u, v1, v2);
        double norm = 1.0/v1S(u);
        uEa1Tv1(u,norm,u);
    }
    
    //random point in square centered on origin with each edge 2*rmax
    public static void randomVector(double[] u, double rmax, Random rand) {
        u[0] = (2.0*rand.nextDouble()-1.0)*rmax;
        u[1] = (2.0*rand.nextDouble()-1.0)*rmax;
    }
    
    public final void randomVector(double[] u, Random rand) {  //random point in entire space
        u[0] = rand.nextDouble()*dimensions[0];
        u[1] = rand.nextDouble()*dimensions[1];
    }
    
}