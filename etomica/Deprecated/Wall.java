package simulate;
import java.io.*;
import java.awt.*;//for Graphics
import java.util.*; //for neighborList
import java.beans.Beans;

public class Wall extends Atom {

    int thickness;
    private boolean vertical, horizontal;
    double[] d = new double[Space.D];   //width and height of wall

    public Wall(Molecule parent, int index) {
        super(parent, index);
        this.setMass(Double.MAX_VALUE);
        setThickness(4);
        setHorizontal(true);
    }
    
    public boolean isVertical() {return vertical;}
    public void setVertical(boolean b) {
        vertical = b;
        if(vertical) {setHorizontal(false);}
    }
    
    public boolean isHorizontal() {return horizontal;}
    public void setHorizontal(boolean b) {
        horizontal = b;
        if(horizontal) {setVertical(false);}
    }
    
    public int getThickness() {return thickness;}
    public void setThickness(int thickness) {this.thickness = thickness;}
 
    public void draw(Graphics g, int[] origin, double scale){
        double toPixels = scale*Phase.TO_PIXELS;
        g.setColor(color);
        int xP = origin[0] + (int)(toPixels*r[0]);
        int yP = origin[1] + (int)(toPixels*r[1]);
        int wP = vertical ? thickness : (int)(toPixels*d[0]);
        int hP = horizontal ? thickness : (int)(toPixels*d[1]);
        g.fillRect(xP,yP,wP,hP);
    }
        
  public double getKineticEnergy() {return 0.0;}
}
