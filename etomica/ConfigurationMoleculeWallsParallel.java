package simulate;

import java.awt.Rectangle;

/**
 * Places all walls parallel to each other, equally spaced
 * 
 * 2-D explicit
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
            x = (double)rect.x/DisplayConfiguration.SIM2PIXELS;
            width = rect.width;
        }
        if(longWall && vertical) {
            y = 0.0;
            height = Integer.MAX_VALUE;
        }
        else {
            y = (double)rect.y/DisplayConfiguration.SIM2PIXELS;
            height = rect.height;
        }
//        int  width = (longWall && phase != null) ? phase.getBounds().width  : rect.width;    //size to phase if longWall, to species otherwise
//        int height = (longWall && phase != null) ? phase.getBounds().height : rect.height;
//        double x = (horizontal && longWall) ? 0.0 : (double)rect.x/Phase.TO_PIXELS;  //put against left or top wall if spanning volume
//        double y = (vertical && longWall)   ? 0.0 : (double)rect.y/Phase.TO_PIXELS;
        double w = (double)width/DisplayConfiguration.SIM2PIXELS;
        double h = (double)height/DisplayConfiguration.SIM2PIXELS;
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
        }                    //2D explicit
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.nextAtom()) {  //equally space all "wall atoms"
            Space2DCell.Vector r = (Space2DCell.Vector)a.coordinate.position();
            if(i==0) {r.x = xyNext;}
            else {r.y = xyNext;}
            xyNext += delta;
            ((AtomType.Wall)a.type).setLength(wh);   //length of wall
            ((AtomType.Wall)a.type).setAngle(angle);
            ((AtomType.Wall)a.type).setTemperature(temperature);
        }
    }
    
    protected void computeDimensions() {
        Molecule m = parentSpecies.makeMolecule();
        initializeCoordinates(m);
        if(horizontal) {
            dim[0] = ((AtomType.Wall)m.firstAtom().type).getLength();
            dim[1] = 0.0;
        }
        else if(vertical) {
            dim[0] = 0.0;
            dim[1] = ((AtomType.Wall)m.firstAtom().type).getLength();
        }
        else {
            //does not handle walls that are neither horizontal nor vertical
        }
    }
}
