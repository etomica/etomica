package etomica.space2d;

import etomica.Default;
import etomica.Phase;
import etomica.Simulation;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
/**
 * Class for implementing slit periodic boundary conditions
 */
public class BoundaryHSlit extends Boundary implements etomica.space.Boundary.Periodic {
    public BoundaryHSlit() {this(Default.BOX_SIZE,Default.BOX_SIZE);}
    public BoundaryHSlit(Phase p) {this(p,Default.BOX_SIZE,Default.BOX_SIZE);}
    public BoundaryHSlit(Phase p, double lx, double ly) {super(p); dimensions.x = lx; dimensions.y = ly; updateDimensions();}
    public BoundaryHSlit(double lx, double ly) {dimensions.x = lx; dimensions.y = ly; updateDimensions();}
    public etomica.space.Boundary.Type type() {return Boundary.HSLIT;}
    private final Vector temp = new Vector();
    private final Vector dimensions = new Vector();
    private final Vector dimensionsCopy = new Vector();
    private final Vector dimensionsHalf = new Vector();
    public final etomica.space.Vector dimensions() {return dimensionsCopy;}
    private final void updateDimensions() {
        dimensionsHalf.Ea1Tv1(0.5,dimensions);
        dimensionsCopy.E(dimensions);
    }
    public etomica.space.Vector randomPosition() {
        temp.x = dimensions.x*Simulation.random.nextDouble(); 
        temp.y = dimensions.y*Simulation.random.nextDouble(); 
        return temp;
    }
    public void nearestImage(etomica.space.Vector dr) {nearestImage((Vector)dr);}
    public void nearestImage(Vector dr) {
       // dr.x -= dimensions.x * ((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x+0.5) : Math.ceil(dr.x/dimensions.x-0.5));
       // dr.y -= dimensions.y * ((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y+0.5) : Math.ceil(dr.y/dimensions.y-0.5));
       // final double dimxHalf = 0.5*dimensions.x;
       // final double dimyHalf = 0.5*dimensions.y;
        while(dr.y > +dimensionsHalf.y) {dr.y -= dimensions.y;}//System.out.println("Space2D here0 3");}
        while(dr.y < -dimensionsHalf.y) {dr.y += dimensions.y;}//System.out.println("Space2D here0 4");}       
    }
    public boolean centralImage(Coordinate c) {return centralImage(c.r);}
    public boolean centralImage(etomica.space.Vector r) {return centralImage((Vector)r);}
    public boolean centralImage(Vector r) {
    /*    while(r.x > dimensions.x) r.x -= dimensions.x;
        while(r.x < 0.0)          r.x += dimensions.x;
        while(r.y > dimensions.y) r.y -= dimensions.y;
        while(r.y < 0.0)          r.y += dimensions.y;
    */  boolean changed = false;
    
        while(r.y > dimensions.y) {r.y -= dimensions.y; changed = true;}//System.out.println("Space2D here 3");}
        while(r.y < 0.0)          {r.y += dimensions.y; changed = true;}//System.out.println("Space2D here 4");}
        return changed;  
        //r.x -= dimensions.x * ((r.x >= 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
        //r.y -= dimensions.y * ((r.y >= 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y-1.0));
    }
    public void inflate(double scale) {
        dimensions.TE(scale); 
        updateDimensions();
        phase().boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
    }
    public void inflate(etomica.space.Vector scale) {
        dimensions.TE(scale); 
        updateDimensions();
        phase().boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
    }
    public void setDimensions(etomica.space.Vector v) {dimensions.E(v); updateDimensions();}
    public double volume() {return dimensions.x * dimensions.y;}
    
    /** 
     * Computes origins for periodic images
     */
    public double[][] imageOrigins(int nShells) {
        int nImages = (2*nShells+1)*(2*nShells+1)-1;
        double[][] origins = new double[nImages][2];
        int k = 0;
        for(int i=-nShells; i<=nShells; i++) {
            for(int j=-nShells; j<=nShells; j++) {
                if(i==0 && j==0) {continue;}
                origins[k][0] = i*dimensions.x;
                origins[k][1] = j*dimensions.y;
                k++;
            }
        }
        return origins;
    }

    /** Returns coordinate shifts needed to draw all images that overflow into central image
     * 0, 1, or 3 shifts may be returned
     */
    int shiftX, shiftY;
    Vector r;
    public float[][] getOverflowShifts(etomica.space.Vector rr, double distance) {
        shiftX = 0; shiftY = 0;
        r = (Vector)rr;
        
        if(r.x-distance < 0.0) {shiftX = +1;}
        else if(r.x+distance > dimensions.x) {shiftX = -1;}
        
        if(r.y-distance < 0.0) {shiftY = +1;}
        else if(r.y+distance > dimensions.y) {shiftY = -1;}
        
        if((shiftX == 0) && (shiftY == 0)) {
          shift = shift0;
        } else if((shiftX != 0) && (shiftY == 0)) {
          shift = new float[1][2];
          shift[0][0] = (float)(shiftX*dimensions.x);
        } else if((shiftX == 0) && (shiftY != 0)) {
          shift = new float[1][2];
          shift[0][1] = (float)(shiftY*dimensions.y);
        } else if((shiftX != 0) && (shiftY != 0)) {
          shift = new float[3][2];
          shift[0][0] = (float)(shiftX*dimensions.x);
          shift[1][1] = (float)(shiftY*dimensions.y);
          shift[2][0] = shift[0][0];
          shift[2][1] = shift[1][1];
        }
        return(shift);
    } //end of getOverflowShifts
}
