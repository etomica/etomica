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
    private double placement;
    
    public ConfigurationMoleculeWallsParallel(){
        setAngle(0);
        setLongWall(false);
        setPlacement(0.0);
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
     * Placement of first wall, as a fraction of the distance from the origin to the end of the phase
     */
    public final double getPlacement() {return placement;}
    public final void setPlacement(double p) {placement = p;}

  /**
   * Sets wall coordinates based on pixel position as obtained by the species' getBounds.
   * Values of x and y coordinates and wall length, are affected by current values of angle and longWall,
   * but these do not in turn alter the original values of the species Bounds
   */
    public void initializeCoordinates(Molecule m) {  //doesn't handle wall that is not either horizontal or vertical
        Rectangle rect = parentSpecies.getBounds();
        Space.Vector d = m.parentPhase().dimensions();
 //       Phase phase = parentSpecies.getParentPhase();
        double x, y;
        double h = d.component(1);
        double w = d.component(0);
        if(longWall) {
            if(horizontal) {
                x = 0.0;
                y = placement*d.component(1);
                w = Double.MAX_VALUE;
            }
            else {//vertical
                x = placement*d.component(0);
                y = 0.0;
                h = Double.MAX_VALUE;
            }
        }
        else {  //finite wall
            x = getBounds().x/DisplayConfiguration.SIM2PIXELS;
            y = getBounds().y/DisplayConfiguration.SIM2PIXELS;
            if(horizontal) {
                w = d.component(0) - x;
            }
            else {
                h = d.component(1) - y;
            }
        }
        int i = 0;
        double delta;
        double xyNext;
        double wh;
        if(horizontal) {
            delta = (h-y)/(m.atomCount-1);
            i = 1;
            xyNext = y;
            wh = w;
        }
        else { //vertical
            delta = (w-x)/(m.atomCount-1);
            i = 0;
            xyNext = x;
            wh = h;
        }                    //2D explicit
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.nextAtom()) {  //equally space all "wall atoms"
            Space.Vector r = a.coordinate.position();
            a.coordinate.momentum().E(0.0);
            r.setComponent(i,xyNext);
            xyNext += delta;
            ((AtomType.Wall)a.type).setLength(wh);   //length of wall
            ((AtomType.Wall)a.type).setAngle(angle);
            ((AtomType.Wall)a.type).setTemperature(temperature);
        }
    }
    
    protected void computeDimensions() {
        if(parentSpecies()==null) return;
        Molecule m = parentSpecies.getMolecule();
//        initializeCoordinates(m);
        if(horizontal) {
//            dim[0] = ((AtomType.Wall)m.firstAtom().type).getLength();
            dim[0] = Double.MAX_VALUE;
            dim[1] = 0.0;
        }
        else if(vertical) {
            dim[0] = 0.0;
            dim[1] = Double.MAX_VALUE;
//            dim[1] = ((AtomType.Wall)m.firstAtom().type).getLength();
        }
        else {
            //does not handle walls that are neither horizontal nor vertical
        }
    }
}
