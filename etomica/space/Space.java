package etomica.space.continuum;
import etomica.*;
import etomica.space.*;
import etomica.space.boundary.*;

public class Space extends etomica.Space implements EtomicaElement {

    public static String version() {return "Space:01.07.11/"+Space.VERSION;}
    public final int D() {return D;}
    
    public Space() {
        this(3);
    }
    public Space(int D) {
        super(D);
    }
    public Space.Vector origin() {return null;}
//    public Space.Vector makeVector() {return new etomica.space.continuum.Vector(D);}
    public Space.Vector makeVector() {
        switch(D) {
            case 3: return new etomica.space.continuum.Vector3D();
            default: return new etomica.space.continuum.Vector(D);
        }
    }
    public Space.Orientation makeOrientation() {return new etomica.space.continuum.Orientation();}
    public Space.Tensor makeTensor() {return new etomica.space.continuum.Tensor();}
    public Space.Tensor makeRotationTensor() {return new etomica.space.continuum.RotationTensor();}
    public Space.Coordinate makeCoordinate(Atom a) {
        if(a.node instanceof AtomTreeNodeGroup) return new etomica.space.CoordinateGroup(this,a);
        else if(a.type instanceof AtomType.Rotator) return new etomica.space.OrientedCoordinate(this,a);
        else return new etomica.space.Coordinate(this,a);
    }
    public Space.CoordinatePair makeCoordinatePair() {return new etomica.space.CoordinatePair(this);}

    //fix this
    public Space.Boundary.Type[] boundaryTypes() {return null;}
//    public Space.Boundary makeBoundary() {return makeBoundary(Boundary.PERIODIC_SQUARE);}  //default
    public Space.Boundary makeBoundary() {return makeBoundary(null);}  //default
    public Space.Boundary makeBoundary(Space.Boundary.Type t) {
/*        if(t == Boundary.NONE) {return new BoundaryNone();}
        else if(t == Boundary.PERIODIC_SQUARE) {return new BoundaryPeriodicSquare();}
   //     else if(t == Boundary.SLIDING_BRICK) return new BoundarySlidingBrick();
        else return null;*/
        return new BoundaryPeriodicCubic(this);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Continuum space");
        return info;
    }

}//end of Space