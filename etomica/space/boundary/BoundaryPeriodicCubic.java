package etomica.space.boundary;
import etomica.*;
import etomica.space.*;

public final class BoundaryPeriodicCubic extends Space.Boundary implements Space.Boundary.Periodic  {
    
/*    public PeriodicCubic() {this(Default.BOX_SIZE,Default.BOX_SIZE,Default.BOX_SIZE);}
    public PeriodicCubic(Phase p) {this(p,Default.BOX_SIZE,Default.BOX_SIZE,Default.BOX_SIZE);}
    public PeriodicCubic(Phase p, double lx, double ly, double lz) {super(p);dimensions.x=lx; dimensions.y=ly; dimensions.z=lz; updateDimensions();}
    public PeriodicCubic(double lx, double ly, double lz) {super();dimensions.x=lx; dimensions.y=ly; dimensions.z=lz; updateDimensions();}
*/

    public final int D;
    
    public BoundaryPeriodicCubic(Space space) {
 /*       temp = (etomica.space.continuum.Vector)space.makeVector();
        dimensions = (etomica.space.continuum.Vector)space.makeVector();
        dimensions.E(Default.BOX_SIZE);
        dimensionsCopy = (etomica.space.continuum.Vector)space.makeVector();
        dimensionsHalf = (etomica.space.continuum.Vector)space.makeVector();
*/        temp = space.makeVector();
        dimensions = space.makeVector();
        dimensions.E(Default.BOX_SIZE);
        dimensionsCopy = space.makeVector();
        dimensionsHalf = space.makeVector();
        D = space.D();
        updateDimensions();
    }

    //fix this
    public Space.Boundary.Type type() {return Space3D.Boundary.PERIODIC_SQUARE;}

    private final Space.Vector temp;
    private final Space.Vector dimensions;
    private final Space.Vector dimensionsCopy;
    private final Space.Vector dimensionsHalf;
//    private final etomica.space.continuum.Vector temp, dimensions, dimensionsCopy, dimensionsHalf;
    public final Space.Vector dimensions() {return dimensionsCopy;}
    
    public Space.Vector randomPosition() {
        temp.setRandomCube();
        temp.PE(0.5);
        temp.TE(dimensions);
        return temp;
    }
    private final void updateDimensions() {
        dimensionsHalf.Ea1Tv1(0.5,dimensions);
        dimensionsCopy.E(dimensions);
    }
    public void nearestImage(Space.Vector dr) {
    /*    while(dr.component(0) > +dimensionsHalf.component(0)) dr.toArray()[0] -= dimensions.component(0);
        while(dr.component(0) < -dimensionsHalf.component(0)) dr.toArray()[0] += dimensions.component(0);
        while(dr.component(1) > +dimensionsHalf.component(1)) dr.toArray()[1] -= dimensions.component(1);
        while(dr.component(1) < -dimensionsHalf.component(1)) dr.toArray()[1] += dimensions.component(1);
        while(dr.component(2) > +dimensionsHalf.component(2)) dr.toArray()[2] -= dimensions.component(2);
        while(dr.component(2) < -dimensionsHalf.component(2)) dr.toArray()[2] += dimensions.component(2);
     */  
        dr.PE(dimensionsHalf);
        dr.mod(dimensions);
        dr.ME(dimensionsHalf);
        /*
        while(dr.x > +dimensionsHalf.x) dr.x -= dimensions.x;
        while(dr.x < -dimensionsHalf.x) dr.x += dimensions.x;
        while(dr.y > +dimensionsHalf.y) dr.y -= dimensions.y;
        while(dr.y < -dimensionsHalf.y) dr.y += dimensions.y;
        while(dr.z > +dimensionsHalf.z) dr.z -= dimensions.z;
        while(dr.z < -dimensionsHalf.z) dr.z += dimensions.z;
        */
    }
    public void centralImage(Coordinate c) {centralImage(c.r);}
    public void centralImage(Space.Vector r) {
        r.mod(dimensions);
     /*   while(r.x > dimensions.x) r.x -= dimensions.x;
        while(r.x < 0.0)          r.x += dimensions.x;
        while(r.y > dimensions.y) r.y -= dimensions.y;
        while(r.y < 0.0)          r.y += dimensions.y;
        while(r.z > dimensions.y) r.z -= dimensions.z;
        while(r.z < 0.0)          r.z += dimensions.z;*/
        //   r.x -= dimensions.x* ((r.x>0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x - 1.0));
        //   r.y -= dimensions.y *((r.y>0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y - 1.0));
        //   r.z -= dimensions.z *((r.z>0) ? Math.floor(r.z/dimensions.z) : Math.ceil(r.z/dimensions.z - 1.0));
    }
    public void inflate(double scale) {dimensions.TE(scale); updateDimensions();}
    public void inflate(Space.Vector scale) {dimensions.TE(scale); updateDimensions();}
    public void setDimensions(Space.Vector v) {dimensions.E(v); updateDimensions();}
    public double volume() {//find a better way
        double[] d = dimensions.toArray();
        double prod = 1.0;
        for(int i=d.length-1; i>=0; i--) {prod *= d[i];}
        return prod;
    }
                
    /**
        * imageOrigins and getOverFlowShifts are both probably incorrect, if they are
        * even completed.  They should definitely be checked before being implemented.
        */
        
        
        //needs to be fixed for 2D/1D
    int shellFormula, nImages, i, j, k, m;
    double[][] origins;
    public double[][] imageOrigins(int nShells) {
        shellFormula = (2 * nShells) + 1;
        nImages = shellFormula*shellFormula*shellFormula-1;
        origins = new double[nImages][D];
        for (k=0,i=-nShells; i<=nShells; i++) {
            for (j=-nShells; j<=nShells; j++) {
                for (m=-nShells; m<=nShells; m++) {
                    if ((i==0 && j==0) && m==0 ) {continue;}
                    origins[k][0] = i*dimensions.component(0);
                    origins[k][1] = j*dimensions.component(1);
                    origins[k][2] = m*dimensions.component(2);
                    k++;
                }
            }
        }
        return origins;
    }
        
        
    //getOverflowShifts ends up being called by the display routines quite often
    //so, in the interest of speed, i moved these outside of the function;
    int shiftX, shiftY, shiftZ;
    Space.Vector r;
    public float[][] getOverflowShifts(Space.Vector rr, double distance) {
        /*
        shiftX = 0; shiftY = 0; shiftZ = 0;
        r = (Vector)rr;
            
        if(r.x-distance < 0.0) {shiftX = +1;}
        else if(r.x+distance > dimensions.component(0)) {shiftX = -1;}
            
        if(r.y-distance < 0.0) {shiftY = +1;}
        else if(r.y+distance > dimensions.component(1)) {shiftY = -1;}
            
        if(r.z-distance < 0.0) {shiftZ = +1;}
        else if(r.z+distance > dimensions.component(2)) {shiftZ = -1;}
              
        if((shiftX == 0) && (shiftY == 0) && (shiftZ == 0)) {
            shift = shift0;
        } else if((shiftX != 0) && (shiftY == 0) && (shiftZ == 0)) {
            shift = new float[1][D];
            shift[0][0] = (float)(shiftX*dimensions.component(0));
        } else if((shiftX == 0) && (shiftY != 0) && (shiftZ == 0)) {
            shift = new float[1][D];
            shift[0][1] = (float)(shiftY*dimensions.component(1));
        } else if((shiftX == 0) && (shiftY == 0) && (shiftZ != 0)) {
            shift = new float[1][D];
            shift[0][2] = (float)(shiftZ*dimensions.component(2));
        } else if((shiftX != 0) && (shiftY != 0) && (shiftZ == 0)) {
            shift = new float[3][D];
            shift[0][0] = (float)(shiftX*dimensions.component(0));
            shift[1][1] = (float)(shiftY*dimensions.component(1));
            shift[2][0] = shift[0][0];
            shift[2][1] = shift[1][1];
        } else if((shiftX != 0) && (shiftY == 0) && (shiftZ != 0)) {
            shift = new float[3][D];
            shift[0][0] = (float)(shiftX*dimensions.component(0));
            shift[1][2] = (float)(shiftZ*dimensions.component(2));
            shift[2][0] = shift[0][0];
            shift[2][2] = shift[1][2];
        } else if((shiftX == 0) && (shiftY != 0) && (shiftZ != 0)) {
            shift = new float[3][D];
            shift[0][1] = (float)(shiftY*dimensions.component(1));
            shift[1][2] = (float)(shiftZ*dimensions.component(2));
            shift[2][1] = shift[0][1];
            shift[2][2] = shift[1][2];
        } else if((shiftX != 0) && (shiftY != 0) && (shiftZ != 0)) {
            shift = new float[7][D];
            shift[0][0] = (float)(shiftX*dimensions.component(0));
            shift[1][1] = (float)(shiftY*dimensions.component(1));
            shift[2][2] = (float)(shiftZ*dimensions.component(2));
            shift[3][0] = shift[0][0];
            shift[3][1] = shift[1][1];
            shift[4][1] = shift[1][1];
            shift[4][2] = shift[2][2];
            shift[5][0] = shift[0][0];
            shift[5][2] = shift[2][2];
            shift[6][0] = shift[0][0];
            shift[6][1] = shift[1][1];
            shift[6][2] = shift[2][2];
        }
            */
        return(shift);
    }//end of getOverflowShifts
}//end of PeriodicCubic
