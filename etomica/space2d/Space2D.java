package simulate.space2D;
//import simulate.Space;
import java.awt.Graphics;
import java.awt.Color;
import java.util.Random;
import simulate.units.*;

public class Space2D extends simulate.Space {
    
    public static final int D = 2;
    public final int D() {return D;}
    
    public double sphereVolume(double r) {return Math.PI*r*r;}  //volume of a sphere of radius r
    public double sphereArea(double r) {return 2.0*Math.PI*r;}  //surface area of sphere of radius r (used for differential shell volume)
    public String[] boundaryTags() {return Boundary.TAGS;}//String descriptions of the boundary options for use in a simulation construction kit
    public simulate.Space.Vector makeVector() {return new Vector();}
    public simulate.Space.Tensor makeTensor() {return new Tensor();}
    public simulate.Space.Coordinate makeCoordinate(Space.Occupant o) {return new Coordinate(o);}
    public simulate.Space.CoordinatePair makeCoordinatePair(Space.Boundary boundary) {return new CoordinatePair((Boundary)boundary);}
    public simulate.Space.Boundary makeBoundary(int b, Phase p) {
        switch(b) {
            case(Boundary.NONE):            return new BoundaryNone(p);
            case(Boundary.PERIODIC_SQUARE): return new BoundaryPeriodicSquare(p);
            case(Boundary.HARD):            return new BoundaryHard(p);
            case(Boundary.SLIDING_BRICK):   return new BoundarySlidingBrick(p);
            default:                        return null;
        }
    }
    public static final double r2(Vector u1, Vector u2, Boundary b) {
        Vector.WORK.x = u1.x - u2.x;
        Vector.WORK.y = u1.y - u2.y;
        b.nearestImage(Vector.WORK);
        return Vector.WORK.x*Vector.WORK.x + Vector.WORK.y*Vector.WORK.y;
    }
}