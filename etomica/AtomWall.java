package simulate;
import java.io.*;
import java.awt.*;//for Graphics
import java.util.*; //for neighborList
import java.beans.Beans;

public class AtomWall extends Atom {

    int thickness;
    private boolean vertical, horizontal;
    double[] r1 = new double[Space.D];  //other end of line in simulation units (first end is denoted r)
    private int theta;   //orientation (degrees) with respect to positive x-axis
    private double length;  //length of line in simulation units
    double pAccumulator;  //accumulated momentum perpendicular to wall, for calculation of pressure

    public AtomWall(Molecule parent, int index) {
        super(parent, index);
        this.setRm(1.0);
        setThickness(4);
        setAngle(0);   //default is horizontal
        setStationary(true);
    }

    public void setAngle(int t) {
        theta = (t <= 360) ? t : (t % 360);
        horizontal = (theta == 0) || (Math.abs(theta) == 180);
        vertical = (Math.abs(theta) == 90) || (Math.abs(theta) == 270);
        computeR1();
    }
    public int getAngle() {return theta;}
        
    public boolean isVertical() {return vertical;}
    public boolean isHorizontal() {return horizontal;}
    
    public int getThickness() {return thickness;}
    public void setThickness(int thickness) {this.thickness = thickness;}
    
    public double getLength() {return length;}
    public void setLength(double length) {
        this.length = length;
        computeR1();
    }
 
    private void computeR1() {
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
            int wP = vertical ? thickness : (int)(toPixels*diameter);
            int hP = horizontal ? thickness : (int)(toPixels*diameter);
            g.fillRect(xP,yP,wP,hP);
        }
    }
        
  public double kineticEnergy() {return 0.0;}
}
