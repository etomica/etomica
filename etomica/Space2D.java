package etomica;

import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomType;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.space.Boundary;
import etomica.space.Coordinate;
import etomica.space.CoordinatePair;
import etomica.space.Orientation;
import etomica.space.RotationTensor;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space2d.BoundaryHSlit;
import etomica.space2d.BoundaryNone;
import etomica.space2d.BoundaryPeriodicSquare;
import etomica.space2d.CoordinateGroup;
import etomica.space2d.OrientedCoordinate;
import etomica.statmech.MaxwellBoltzmann;

//centralImage not updated to molecule form as in Space3D

 /* History of changes
  * 09/01/02 (DAK) added accelerateTo method to Coordinate
  *                changed CoordinateGroup.randomizeMomentum to not enforce zero COM momentum
  * 09/05/02 (DAK) fixed error in accelerateTo (still probably does not do what one expects
  *                if accelerating to nonzero momentum).
  * 01/04/03 (SKK/DAK) added HSlit Boundary class for horizontal-slit PBC.
  * 01/12/03 (JKS/DAK) corrected error in Vector.transform, where updated xyz
  * 07/10/03 (DAK) added resetV method to CoordinatePair
  * 08/27/03 (DAK) added isZero method to Vector
  * 08/29/03 (DAK) implemented centralImage(Space.Coordinate) in Boundary
  * 12/09/03 (DAK) added setRandomInSphere method in Vector
  * 01/22/04 (DAK) added positionCOM and translateCOMTo to CoordinateGroup;
  * redefined position() in CoordinateGroup to be first-atom position (as it has
  * been for Space3D for some time now).
  */
public class Space2D extends Space implements EtomicaElement {
    
    public static final int D = 2;
    public final int D() {return D;}
    public final int powerD(int n) {return n*n;}
    public final double powerD(double a) {return a*a;}
    public static final Vector ORIGIN = new Vector();
    public final Vector origin() {return ORIGIN;}
    public static final Space2D INSTANCE = new Space2D();
    
    public Space2D() {super(2);}
    
    public double sphereVolume(double r) {return Math.PI*r*r;}  //volume of a sphere of radius r
    public double sphereArea(double r) {return 2.0*Math.PI*r;}  //surface area of sphere of radius r (used for differential shell volume)
    public Vector makeVector() {return new Vector();}
    public Orientation makeOrientation() {return new Orientation();}
    public Tensor makeTensor() {return new Tensor();}
    public Tensor makeRotationTensor() {return new RotationTensor();}
    public Coordinate makeCoordinate(Atom a) {
        if(a.node instanceof AtomTreeNodeGroup) return new CoordinateGroup(a);
        else if(a.type instanceof AtomType.Rotator) return new OrientedCoordinate(a);
        else {return new Coordinate(a);}
    }
    public CoordinatePair makeCoordinatePair() {return new CoordinatePair();}
    
    public Boundary.Type[] boundaryTypes() {return Boundary.TYPES;}
    public Boundary makeBoundary() {return makeBoundary(Boundary.PERIODIC_SQUARE);}  //default
    public Boundary makeBoundary(Boundary.Type t) {
        if(t == Boundary.NONE) {return new BoundaryNone();}
        else if(t == Boundary.PERIODIC_SQUARE) {return new BoundaryPeriodicSquare();}
        else if(t == Boundary.HSLIT) {return new BoundaryHSlit();}
 //       else if(t == Boundary.SLIDING_BRICK) return new BoundarySlidingBrick();
        else return null;
    }
    
    public int[] makeArrayD(int i) {return new int[] {i, i};}
    public double[] makeArrayD(double d) {return new double[] {d, d};}
 
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Two-dimensional space");
        return info;
    }

    public static final double r2(Vector u1, Vector u2, Boundary b) {
        Vector.WORK.x = u1.x - u2.x;
        Vector.WORK.y = u1.y - u2.y;
        b.nearestImage(Vector.WORK);
        return Vector.WORK.x*Vector.WORK.x + Vector.WORK.y*Vector.WORK.y;
    }
  
  /*  public static final class BoundarySlidingBrick extends BoundaryPeriodicSquare {
        private double gamma = 0.0;
        private double delvx;
        private IntegratorMD.ChronoMeter timer;
        public BoundarySlidingBrick() {super();}
        public BoundarySlidingBrick(Phase p) {super(p);}
        public Space.Boundary.Type type() {return Boundary.SLIDING_BRICK;}
        public void setShearRate(double g) {gamma = g; computeDelvx();}
        public double getShearRate() {return gamma;}
        private void computeDelvx() {delvx = gamma*dimensions.y;}
        
        public void setTimer(IntegratorMD.ChronoMeter t) {timer = t;}
        
        public void nearestImage(Vector dr) {
            double delrx = delvx*timer.currentValue();
            double cory = ((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y+0.5) : Math.ceil(dr.y/dimensions.y-0.5));
            dr.x -= cory*delrx;
            dr.x -= dimensions.x * ((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x+0.5) : Math.ceil(dr.x/dimensions.x-0.5));
            dr.y -= dimensions.y * cory;
        }
        public void centralImage(Vector r) {
            double delrx = delvx*timer.currentValue();
            double cory = ((r.y >= 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y-1.0));
//            if(cory != 0.0) System.out.println(delrx*cory);
            r.x -= cory*delrx;
            r.x -= dimensions.x * ((r.x >= 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
            r.y -= dimensions.y * cory;
        }
        public void centralImage(Coordinate c) {
            Vector r = c.r;
            double cory = ((r.y > 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y-1.0));
            double corx = ((r.x > 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
            if(corx==0.0 && cory==0.0) return;
            double delrx = delvx*timer.currentValue();
            Vector p = c.p;
            r.x -= cory*delrx;
            r.x -= dimensions.x * corx; 
            r.y -= dimensions.y * cory;
            p.x -= cory*delvx;
        }
        
        public double[][] imageOrigins(int nShells) {
            int nImages = (2*nShells+1)*(2*nShells+1)-1;
            double[][] origins = new double[nImages][D];
            int k = 0;
            for(int i=-nShells; i<=nShells; i++) {
                for(int j=-nShells; j<=nShells; j++) {
                    if(i==0 && j==0) {continue;}
                    origins[k][0] = i*dimensions.x + j*delvx*timer.currentValue();
                    origins[k][1] = j*dimensions.y;
                    k++;
                }
            }
            return origins;
        }//end of imageOrigins
    }//end of BoundarySlidingBrick */
            
}//end of Space2D
