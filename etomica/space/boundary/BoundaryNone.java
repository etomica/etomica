package etomica.space.boundary;
import etomica.Space;
import etomica.Default;

public final class BoundaryNone extends Space.Boundary {
    
    private final Space.Vector temp;
    public final Space.Vector dimensions;
    private final int D;
    //fix this
//    public static final Space.Boundary.Type NONE = new Space.Boundary.Type("None");
    public static final Space.Boundary.Type NONE = etomica.Space3D.Boundary.NONE;
    
    public BoundaryNone(Space space) {
        super();
        temp = space.makeVector();
        dimensions = space.makeVector();
        dimensions.E(Default.BOX_SIZE);
        D = space.D();
    }
    
    public Space.Boundary.Type type() {return NONE;}
    public final Space.Vector dimensions() {return dimensions;}
    public void nearestImage(Space.Vector dr) {}
    public boolean centralImage(Space.Vector r) {return false;}
    public boolean centralImage(Space.Coordinate c) {return false;}
    public double volume() {//find a better way
        double[] d = dimensions.toArray();
        double prod = 1.0;
        for(int i=d.length-1; i>=0; i--) {prod *= d[i];}
        return prod;
    }
    public void inflate(double s) {dimensions.TE(s);}
    public void inflate(Space.Vector s) {dimensions.TE(s);}
    public void setDimensions(Space.Vector v) {dimensions.E(v);}
    public double[][] imageOrigins(int nShells) {return new double[0][D];}
    public float[][] getOverflowShifts(Space.Vector rr, double distance) {return shift0;}
    public Space.Vector randomPosition() {
        temp.setRandomCube();
        temp.PE(0.5);
        temp.TE(dimensions);
        return temp;
    }
}//end of None
