package simulate;

import java.awt.Rectangle;

/**
 * Places all walls parallel to each other, equally spaced
 */

public class ConfigurationMoleculeWallsParallel extends ConfigurationMolecule {
    
    private int angle;
    private boolean horizontal, vertical;
    private boolean longWall;  //If true, specifies that the wall extends the whole length of the simulation volume
    private double temperature = 300.;
    
    public ConfigurationMoleculeWallsParallel(){
        setAngle(0);
        setLongWall(false);
    }
      
    public final int getAngle() {return angle;}
    public final void setAngle(int t) {
        t = (Math.abs(t) > 45) ? 90 : 0;  //For now, allow only values for vertical or horizontal walls
        angle = (t <= 360) ? t : (t % 360);
        horizontal = (angle == 0) || (Math.abs(angle) == 180);
        vertical = (Math.abs(angle) == 90) || (Math.abs(angle) == 270);
        initializeCoordinates();
    }
    
    public final double getTemperature() {return temperature;}
    public final void setTemperature(double t) {
        temperature = t;
        initializeCoordinates();
    }
    
    public final boolean isLongWall() {return longWall;}
    public final void setLongWall(boolean s) {
        longWall = s;
        initializeCoordinates();
    }

  /**
   * Sets wall coordinates based on pixel position as obtained by the species' getBounds.
   * Values of x and y coordinates and wall length, are affected by current values of angle and longWall,
   * but these do not in turn alter the original values of the species Bounds
   */
    public void initializeCoordinates(Molecule m) {  //doesn't handle wall that is not either horizontal or vertical
        Rectangle rect = parentSpecies.getBounds();
 //       Phase phase = parentSpecies.getParentPhase();
        double x, y;
        int width, height;
        if(longWall && horizontal) {
            x = 0.0;
            width = Integer.MAX_VALUE;
        }
        else {
            x = (double)rect.x/Phase.TO_PIXELS;
            width = rect.width;
        }
        if(longWall && vertical) {
            y = 0.0;
            height = Integer.MAX_VALUE;
        }
        else {
            y = (double)rect.y/Phase.TO_PIXELS;
            height = rect.height;
        }
//        int  width = (longWall && phase != null) ? phase.getBounds().width  : rect.width;    //size to phase if longWall, to species otherwise
//        int height = (longWall && phase != null) ? phase.getBounds().height : rect.height;
//        double x = (horizontal && longWall) ? 0.0 : (double)rect.x/Phase.TO_PIXELS;  //put against left or top wall if spanning volume
//        double y = (vertical && longWall)   ? 0.0 : (double)rect.y/Phase.TO_PIXELS;
        double w = (double)width/Phase.TO_PIXELS;
        double h = (double)height/Phase.TO_PIXELS;
        int i = 0;
        double delta = 0.0;
        double xyNext = x;
        double wh = h;
        if(horizontal) {
            delta = (h-y)/(m.nAtoms-1);
            i = 1;
            xyNext = y;
            wh = w;
        }
        else if(vertical) {
            delta = (w-x)/(m.nAtoms-1);
            i = 0;
            xyNext = x;
            wh = h;
        }
        for(AtomC a=(AtomC)m.firstAtom(); a!=m.terminationAtom(); a=(AtomC)a.getNextAtom()) {  //equally space all "wall atoms"
            Space.uEa1(a.r,0.0);
            a.r[i] = xyNext;
            xyNext += delta;
            ((AtomWall)a).setLength(wh);   //length of wall
            ((AtomWall)a).setAngle(angle);
            ((AtomWall)a).setTemperature(temperature);
        }
    }
    
    protected void computeDimensions() {
        if(horizontal) {
            dim[0] = ((AtomHardWall)parentSpecies.firstAtom()).getLength();
            dim[1] = 0.0;
        }
        else if(vertical) {
            dim[0] = 0.0;
            dim[1] = ((AtomHardWall)parentSpecies.firstAtom()).getLength();
        }
        else {
            //does not handle walls that are neither horizontal nor vertical
        }
    }
}
