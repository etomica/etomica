package simulate;
import java.io.*;
import java.awt.*;//for Graphics
import java.util.*; //for neighborList
import java.beans.Beans;

public class AtomHardWall extends AtomHard implements AtomWall {

    int thickness;
    private boolean vertical, horizontal;
    double[] f = new double[Space.D];   //force on wall
    double[] r1 = new double[Space.D];  //other end of line in simulation units (first end is denoted r)
    protected int theta;   //orientation (degrees) with respect to positive x-axis, taking r at the origin
    protected double length;  //length of line in simulation units
    protected double pAccumulator;  //accumulated momentum perpendicular to wall, for calculation of pressure
    protected double qAccumulator;  //accumulated energy absorbed by wall, for calculation of heat transfer
    protected double temperature;
    protected boolean adiabatic;

    /**
     * Constructs an atom with no initialization if parent is null; otherwise constructs atom with default atomIndex = 0.  
     * Expected use of such an Atom is for the construction of other Atoms via makeAtom method
     */
    public AtomHardWall(Molecule parent) {
        this(parent,0);
    }
    
    public AtomHardWall(Molecule parent, int index) {
        super(parent, index);
        if(parent != null) {
            this.setRm(1.0);
            setThickness(4);
            setAngle(0);   //default is horizontal
            setStationary(true);
            setTemperature(300.);
            setAdiabatic(true);
            zeroAccumulators();
        }
    }

    public Atom makeAtom(Molecule m, int i) {return new AtomHardWall(m,i);}

    public void setAngle(int t) {
        theta = (t <= 360) ? t : (t % 360);
        horizontal = (theta == 0) || (Math.abs(theta) == 180);
        vertical = (Math.abs(theta) == 90) || (Math.abs(theta) == 270);
        computeR1();
    }
    public final int getAngle() {return theta;}
        
    public final boolean isVertical() {return vertical;}
    public final boolean isHorizontal() {return horizontal;}
    
    public final int getThickness() {return thickness;}
    public final void setThickness(int thickness) {this.thickness = thickness;}
    
    public final double getLength() {return length;}
    public final void setLength(double length) {
        this.length = length;
        computeR1();
    }
    
    public final void zeroAccumulators() {
        zeroQAccumulator();
        zeroPAccumulator();
    }
    public final void zeroQAccumulator() {qAccumulator = 0.0;}
    public final void zeroPAccumulator() {pAccumulator = 0.0;}
    public final void accumulateP(double value) {pAccumulator += value;}
    public final void accumulateQ(double value) {qAccumulator += value;}
    public final double getAccumulatedP() {return pAccumulator;}
    public final double getAccumulatedQ() {return qAccumulator;}
    
    public final double getTemperature() {return temperature;}
    public final void setTemperature(int t) {setTemperature((double)t);}  //for connection to sliders, etc.
    public final void setTemperature(double t) {temperature = t;}
    
    public final boolean isAdiabatic() {return adiabatic;}
    public final void setAdiabatic(boolean a) {adiabatic = a;}
 
    private final void computeR1() {
        r1[0] = r[0] + length * Math.cos((double)theta);
        r1[1] = r[1] + length * Math.sin((double)theta);
    }
    public void draw(Graphics g, int[] origin, double scale){
        double toPixels = scale*Phase.TO_PIXELS;
        g.setColor(color);
        int xP = origin[0] + (int)(toPixels*r[0]);
        int yP = origin[1] + (int)(toPixels*r[1]);
        if(!horizontal && !vertical) {
            int x1 = xP + (int)(toPixels*r1[0]);
            int y1 = yP + (int)(toPixels*r1[1]);
            g.drawLine(xP, yP, x1, y1);
        }
        else {    
            int wP = vertical ? thickness : (int)(toPixels*length);
            int hP = horizontal ? thickness : (int)(toPixels*length);
            g.fillRect(xP,yP,wP,hP);
        }
    }
        
  public double kineticEnergy() {return 0.0;}
}
