package simulate;

/**
 * Places all walls along the boundary of phase
 * Keeps placing walls around phase, until all are positioned
 * Direction of placement can be specified clockwise or counterclockwise
 */

public class ConfigurationMoleculeWallsBoundary extends ConfigurationMolecule {
    
    private double temperature = 300.;
    private boolean fillClockwise = true;
    private int firstWallPosition = Constants.NORTH;
    
    public ConfigurationMoleculeWallsBoundary(){
    }

    public final double getTemperature() {return temperature;}
    public final void setTemperature(double t) {
        temperature = t;
        initializeCoordinates();
    }
    
    public final int getFirstWallPosition() {return firstWallPosition;}
    public final void setFirstWallPosition(int i) {
        firstWallPosition = i;
        initializeCoordinates();
    }
    
    public boolean isFillClockwise() {return fillClockwise;}
    public void setFillClockwise(boolean b) {
        fillClockwise = b;
        initializeCoordinates();
    }
    
   
    public void initializeCoordinates(Molecule m) {
        if(parentSpecies.getParentPhase() == null) {return;}
        AtomC first = (AtomC)m.firstAtom();
        AtomC last = (AtomC)m.lastAtom();
        switch(firstWallPosition) {  //one method call places first and all subsequent atoms
            case Constants.NORTH: placeNorth(first, last); break;
            case Constants.SOUTH: placeSouth(first, last); break;
            case Constants.EAST:   placeEast(first, last); break;
            case Constants.WEST:   placeWest(first, last); break;
        }
    }
    
    private void placeNorth(AtomC wf, AtomC wl) {
        if(wf == null) {return;}
        Space.uEa1(wf.r,0.0);
        ((AtomWall)wf).setAngle(0);
        if(wf == wl) {return;}
        if(fillClockwise) {placeEast((AtomC)wf.getNextAtom(),wl);}
        else {placeWest((AtomC)wf.getNextAtom(),wl);}
    }
        
    private void placeSouth(AtomC wf, AtomC wl) {
        if(wf == null) {return;}
        Space.uEa1(wf.r,0.0);
        wf.r[1] = (parentSpecies.getParentPhase().getSize().height-((AtomWall)wf).getThickness())/Phase.TO_PIXELS;
        ((AtomWall)wf).setAngle(0);
        if(wf == wl) {return;}
        if(fillClockwise) {placeWest((AtomC)wf.getNextAtom(),wl);}
        else {placeEast((AtomC)wf.getNextAtom(),wl);}
    }
        
    private void placeEast(AtomC wf, AtomC wl) {
        if(wf == null) {return;}
        Space.uEa1(wf.r,0.0);
        wf.r[0] = (parentSpecies.getParentPhase().getSize().width-((AtomWall)wf).getThickness())/Phase.TO_PIXELS;
        ((AtomWall)wf).setAngle(90);
        if(wf == wl) {return;}
        if(fillClockwise) {placeSouth((AtomC)wf.getNextAtom(),wl);}
        else {placeNorth((AtomC)wf.getNextAtom(),wl);}
    }
        
    private void placeWest(AtomC wf, AtomC wl) {
        if(wf == null) {return;}
        Space.uEa1(wf.r,0.0);
        ((AtomWall)wf).setAngle(90);
        if(wf == wl) {return;}
        if(fillClockwise) {placeNorth((AtomC)wf.getNextAtom(),wl);}
        else {placeSouth((AtomC)wf.getNextAtom(),wl);}
    }
        
    protected void computeDimensions() {
        dim.E(1.0);  //no meaningful choice possible here.
    }
}