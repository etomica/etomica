package simulate;
import java.io.*;
import java.awt.*;  //subclasses Component; for Graphics and Color in draw method
import java.beans.*;
import java.util.*;

public class Space extends Component {

    public static final int D = 2;                       //dimension of space
    public boolean periodic;
    protected final double[] dimensions = new double[D]; //size: width, height in simulation units
    protected final int[] drawSize = new int[D];         //size of space, in pixels
    protected double scale = 1.0;
    private static final double[] work = new double[D];       //temporary work vector
    protected final double[][] shift0 = new double[0][D];
    protected Phase parentPhase;
    public double volume;
    
    protected int nShells = 0;
    protected int nImages;
    protected int[][] origins;
        
    public Space() {
        dimensions[0] = 1.0;
        dimensions[1] = 1.0;
        computeVolume();
        periodic = false;
    }
    
    public void setParentPhase(Phase p) {
        parentPhase = p;
        setDimensions(0,dimensions[0]);
    }
    
    public final boolean isPeriodic() {return periodic;}
    
    public void draw(Graphics g, int[] origin) {
        if(isVisible()) {
            g.setColor(Color.gray.brighter());
            g.drawRect(origin[0],origin[1],drawSize[0]-1,drawSize[1]-1);
        }
    }
    
    public double getScale() {return scale;}
    public void setScale(double s, int n) {   //likely to override
        if(s>0 && n>=0) {
          scale = s;
          nShells = n;
          computeDrawSize();
        }
    }
    
    protected void computeDrawSize() {
        drawSize[0] = Phase.toPixels(scale*dimensions[0]);
        drawSize[1] = Phase.toPixels(scale*dimensions[1]);
        resetOrigins(nShells);
    }
    
    public int[] getDrawSize() {return drawSize;}

    public double[] getDimensions() {return dimensions;}
    
    public double getDimensions(int i) {return dimensions[i];}
    public void setDimensions(int i, double dim) {
        dimensions[i] = dim;
        if(parentPhase != null) {parentPhase.resetTO_PIXELS();}
        computeDrawSize();
        computeVolume();
    }
    
    public void inflate(double scale) {
        setDimensions(0,scale*dimensions[0]);
        setDimensions(1,scale*dimensions[1]);
    }
    
    protected void computeVolume() {
        volume = dimensions[0]*dimensions[1];
    }
    
    protected void resetOrigins(int n) {    //likely to override
        nShells = n;
        nImages = 0;
        origins = new int[nImages][];
    }
        
    public int[][] getImageOrigins(int n) {
        if(n != nShells) {resetOrigins(n);}
        return origins;
    }
    
    public double[][] getOverflowShifts(double[] r, double distance) {return shift0;}  //called only if periodic
        
    public void setXDimension(double dim) {setDimensions(0, dim);} //used for property list editor
    public double getXDimension() {return dimensions[0];}
    public void setYDimension(double dim) {setDimensions(1,dim);}
    public double getYDimension() {return dimensions[0];}
    
    public void repositionMolecules(Species species) {return;}  

    public void paint(Graphics g) {
        if(Beans.isDesignTime()) {
            g.setColor(Color.yellow);
            g.fillRect(0,0,getSize().width,getSize().height);
        }
    }
 
    /* The remainder of the methods do vector arithmetic.
    
       In the argument list, "u" represents the "left-hand side" if the 
       method calculation results in a vector; the array passed in this
       position will have its values updated according to the calculation.
       If the calculation results in a scalar, this is the return value
       of the method and no "u" is present in the argument list.
       
       v1, v2, etc. represent vectors appearing on the right-hand side
       of the vector equation, while a1, a2, etc. represent scalars 
       appearing there.
       
       Method names are coded to suggest how the vector equation
       would be written using standard notation.  Thus:
       
       uE --> "u equals"; returns a vector in the first argument
        P --> "plus"
        M --> "minus"
        T --> "times"
        D --> "dot"; inner product
        S --> "squared"; vector square magnitude
       PE --> "plus equals"  (+=); increments u
       ME --> "minus equals" (-=)
       underscores represent parentheses ("left paren" omitted if would be first character of method name)
    */  
    
    public static void uEa1(double[] u, double a1) {
        u[0] = u[1] = a1;
    }
    
    public static void uEv1(double[] u, double[] v1) {
        u[0] = v1[0];
        u[1] = v1[1];
    }
    
    public static double v1Dv2(double[] v1, double[] v2) {
        return (v1[0]*v2[0] + v1[1]*v2[1]);
    }
    
    public static double v1S(double[] v1) {
        return (v1[0]*v1[0] + v1[1]*v1[1]);
    }
    
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
        
    public static void uTEa1(double[] u, double a1) {
        u[0] *= a1;
        u[1] *= a1;
        return;
    }
    
    public static void uTEa1(double[] u, int a1) {
        u[0] *= a1;
        u[1] *= a1;
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