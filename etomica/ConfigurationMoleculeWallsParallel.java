package simulate;

/**
 * Places all atoms of each molecule in a straight line along the x-axis
 */

public class ConfigurationMoleculeWallsParallel extends ConfigurationMolecule {
    
    private double theta;
    private boolean horizontal, vertical;
    private boolean spanVolume;
    
    public ConfigurationMoleculeWallsParallel(){
        setAngle(0);
        setSpanVolume(true);
    }
      
    public final int getAngle() {return theta;}
    public final void setAngle(int t) {
        t = (Math.abs(t) > 45) ? 90 : 0;  //For now, allow only values for vertical or horizontal walls
        theta = (t <= 360) ? t : (t % 360);
        horizontal = (theta == 0) || (Math.abs(theta) == 180);
        vertical = (Math.abs(theta) == 90) || (Math.abs(theta) == 270);
    }
    
    public final boolean getSpanVolume() {return spanVolume;}
    public final void setSpanVolume(boolean s) {spanVolume = s;}

  /**
   * Sets wall coordinates based on pixel position as obtained by the species' getBounds.
   * Values of x and y coordinates and wall length, are affected by current values of angle and spanVolume,
   * but these do not in turn alter the original values of the species Bounds
   */
    public void initializeCoordinates(Molecule m) {
        Rectangle rect = parentSpecies.getBounds();
        double x = (double)rect.x/Phase.TO_PIXELS;
        double y = (double)rect.y/Phase.TO_PIXELS;
        double width = (double)rect.width/Phase.TO_PIXELS;
        double height = (double)rect.height/Phase.TO_PIXELS;
        if(horizontal) {
            double delta = (height-y)/(m.nAtoms-1);
            Space.uEa1(m.firstAtom().r,0.0);
            Atom a = m.firstAtom();
            a.r[1] = y;
            m.firstAtom.setDiameter(width);
            for(Atom a=a.getNextAtom(); a!=m.terminationAtom(); a=a.getNextAtom()) {
                Space.uEa1(a.r,0.0);
                a.r[1] = a.getPreviousAtom().r[1] + delta;
                a.setDiameter(width);
            }
        }
        else if(vertical) {
            double delta = (width-x)/(m.nAtoms-1);
            Space.uEa1(m.firstAtom().r,0.0);
            Atom a = m.firstAtom();
            a.r[0] = x;
            m.firstAtom.setDiameter(height);
            for(Atom a=a.getNextAtom(); a!=m.terminationAtom(); a=a.getNextAtom()) {
                Space.uEa1(a.r,0.0);
                a.r[0] = a.getPreviousAtom().r[0] + delta;
                a.setDiameter(height);
            }
            
            //doesn't handle wall that is not either horizontal or vertical
        }
    }
    
    
    protected void computeDimensions() {
        dim[1] = 0.0;
        Molecule m = parentSpecies.firstMolecule();  //a typical molecule
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.getNextAtom()) {
            dim[1] = Math.max(dim[1], a.getDiameter());  //width is that of largest atom
        }
        dim[0] = 0.5*(m.firstAtom().getDiameter() + m.lastAtom().getDiameter()) + (m.nAtoms-1) * bondLength;
    }
}
