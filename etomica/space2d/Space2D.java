package etomica.space2d;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Space;
import etomica.space.Boundary;

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
    public static final Vector2D ORIGIN = new Vector2D();
    public final etomica.space.Vector origin() {return ORIGIN;}
    public static final Space2D INSTANCE = new Space2D();
    
    public Space2D() {super(2);}
    
    public double sphereVolume(double r) {return Math.PI*r*r;}  //volume of a sphere of radius r
    public double sphereArea(double r) {return 2.0*Math.PI*r;}  //surface area of sphere of radius r (used for differential shell volume)
    public etomica.space.Vector makeVector() {return new Vector2D();}
    public etomica.space.Orientation makeOrientation() {return new Orientation();}
    public etomica.space.Tensor makeTensor() {return new Tensor2D();}
    public etomica.space.Tensor makeRotationTensor() {return new RotationTensor();}
    
    public int[] makeArrayD(int i) {return new int[] {i, i};}
    public double[] makeArrayD(double d) {return new double[] {d, d};}
 
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Two-dimensional space");
        return info;
    }

    public static final double r2(Vector2D u1, Vector2D u2, Boundary b) {
        Vector2D.WORK.x = u1.x - u2.x;
        Vector2D.WORK.y = u1.y - u2.y;
        b.nearestImage(Vector2D.WORK);
        return Vector2D.WORK.x*Vector2D.WORK.x + Vector2D.WORK.y*Vector2D.WORK.y;
    }
            
}//end of Space2D
