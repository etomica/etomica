package simulate;

import java.awt.Rectangle;

/**
 * Places all atoms of each molecule in a straight line along the x-axis
 */

public class ConfigurationMoleculeWallsParallel extends ConfigurationMolecule {
    
    private int angle;
    private boolean horizontal, vertical;
    private boolean spanVolume;
    private double temperature = 300.;
    
    public ConfigurationMoleculeWallsParallel(){
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
    
    public final boolean getSpanVolume() {return spanVolume;}
    public final void setSpanVolume(boolean s) {
        spanVolume = s;
        initializeCoordinates();
    }

  /**
   * Sets wall coordinates based on pixel position as obtained by the species' getBounds.
   * Values of x and y coordinates and wall length, are affected by current values of angle and spanVolume,
   * but these do not in turn alter the original values of the species Bounds
   */
    public void initializeCoordinates(Molecule m) {  //doesn't handle wall that is not either horizontal or vertical
        Rectangle rect = parentSpecies.getBounds();
        Phase phase = parentSpecies.getParentPhase();
        int  width = (spanVolume && phase != null) ? phase.getBounds().width  : rect.width;    //size to phase if spanVolume, to species otherwise
        int height = (spanVolume && phase != null) ? phase.getBounds().height : rect.height;
        double x = (horizontal && spanVolume) ? 0.0 : (double)rect.x/Phase.TO_PIXELS;  //put against left or top wall if spanning volume
        double y = (vertical && spanVolume)   ? 0.0 : (double)rect.y/Phase.TO_PIXELS;
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
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.getNextAtom()) {  //equally space all "wall atoms"
            Space.uEa1(a.r,0.0);
            a.r[i] = xyNext;
            xyNext += delta;
            a.setDiameter(wh);
            ((AtomWall)a).setAngle(angle);
            ((AtomWall)a).setTemperature(temperature);
        }
    }
    
    protected void computeDimensions() {
        if(horizontal) {
            dim[0] = ((AtomWall)parentSpecies.firstAtom()).getLength();
            dim[1] = 0.0;
        }
        else if(vertical) {
            dim[0] = 0.0;
            dim[1] = ((AtomWall)parentSpecies.firstAtom()).getLength();
        }
        else {
            //does not handle walls that are neither horizontal nor vertical
        }
    }
}
