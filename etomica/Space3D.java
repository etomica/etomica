package etomica;

import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomType;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.space3d.*;
import etomica.statmech.MaxwellBoltzmann;

/**
 *
 * @author Rob Riggleman
 * @author David Kofke
 */
 
 /* History of changes
  * 7/24/02 started recording change history
  * 7/24/02 (DW) modified RotationTensor
  *         (DW) change to CoordinatePair to use vector operators, needed to ensure proper working with CoordinateGroup
  *         (?) unknown changes to Orientation
  * 09/01/02 (DAK) added accelerateTo method to Coordinate
  *                changed CoordinateGroup.randomizeMomentum to not enforce zero COM momentum
  * 09/05/02 (DAK) fixed error in accelerateTo (still probably does not do what one expects
  *                if accelerating to nonzero momentum).
  * 01/12/03 (JKS/DAK) corrected error in Vector.transform, where updated xyz
  * values were being used prematurely
  * 01/31/03 (JKS/DAK) modifications to make thread-safe
  * 07/10/03 (DAK) added resetV method to CoordinatePair
  * 08/13/03 (DAK) added massSum check in CoordinateGroup.position()
  * 08/25/03 (DAK/DW) revised CoordinateGroup.position() to return position of
  * first atom, rather than center-of-mass position
  * 08/27/03 (DAK) added isZero and EModShift methods to Vector
  * 08/28/03 (DAK) added mod2, Mv1Pv2Squared, Mv1Squared methods to Vector;
  * added nearestImage(dr, shift) method to Boundary with
  * nontrivial implementation in BoundaryPeriodicCubic
  * 08/29/03 (DAK) implemented centralImage(Space.Coordinate) in Boundary
  * 08/29/03 (DAK) translateBy(Space.Vector) invokes translateBy(Vector)
  * instead of repeating the code it it, then in CoordinateGroup translateBy
  * (Vector) overrides Coordinate.translateBy(Vector) method.
  * 12/09/03 (DAK) added setRandomInSphere method in Vector
  * 01/22/04 (DAK) added positionCOM and translateCOMTo to CoordinateGroup
  */

public class Space3D extends Space implements EtomicaElement {

    public static final int D = 3;
    public final int D() {return D;}
    public final int powerD(int n) {return n*n*n;}
    public final double powerD(double a) {return a*a*a;}
    public int[] makeArrayD(int i) {return new int[] {i, i, i};}
    public double[] makeArrayD(double d) {return new double[] {d, d, d};}
    
    public static final Vector ORIGIN = new Vector();
    public final etomica.space.Vector origin() {return ORIGIN;}
    public static final Space3D INSTANCE = new Space3D();
    
    public Space3D() {super(3);}
    
    public double sphereVolume(double r) {return (Math.PI*4.0*r*r*r/3.0);}
    public double sphereArea(double r)  {return (Math.PI*4*r*r);}
    public etomica.space.Vector makeVector() {return new Vector();}
    public etomica.space.Orientation makeOrientation() {return new Orientation();}
    public etomica.space.Tensor makeTensor() {return new Tensor();}
    public etomica.space.Tensor makeRotationTensor() {return new RotationTensor();}
    public etomica.space.Coordinate makeCoordinate(Atom a) {
        if(a.node instanceof AtomTreeNodeGroup) return new CoordinateGroup(a);
        else if(a.type instanceof AtomType.Rotator) return new OrientedCoordinate(a);
        else return new Coordinate(a);
    }
    public etomica.space.CoordinatePair makeCoordinatePair() {return new CoordinatePair();}

    public etomica.space.Boundary.Type[] boundaryTypes() {return Boundary.TYPES;}
    public etomica.space.Boundary makeBoundary() {return makeBoundary(Boundary.PERIODIC_SQUARE);}  //default
    public etomica.space.Boundary makeBoundary(etomica.space.Boundary.Type t) {
        if(t == Boundary.NONE) {return new BoundaryNone();}
        else if(t == Boundary.PERIODIC_SQUARE) {return new BoundaryPeriodicSquare();}
   //     else if(t == Boundary.SLIDING_BRICK) return new BoundarySlidingBrick();
        else return null;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Three-dimensional space");
        return info;
    }

    public static final double r2(Vector u1, Vector u2, Boundary b) {
    	return r2(u1, u2, b, new Vector());
    }
    public static final double r2(Vector u1, Vector u2, Boundary b, Vector work) {
    	work.Ev1Mv2(u1, u2);
        b.nearestImage(work);
        return work.squared();
    }  
    
    public static void main (String[] args) {
        Vector r1 = new Vector(2,2,3);
        System.out.println("r1_before" + r1.toString());
        Tensor tensor = new Tensor(new double[] {1,2,0,1,1,2,0,0,1});
        RotationTensor tensor2 = new RotationTensor();
        //r1.transform(tensor2);
        tensor2.E(tensor);
        System.out.println("tensor2_before " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
        System.out.println();
        
        r1.transform(tensor2);
        System.out.println("r1_transform(tensor2)" + r1.toString());
        tensor2.invert();
        System.out.println("tensor2_invert " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
        tensor2.setAxial(1, 2*Math.PI);
        System.out.println("tensor2_rotate_360 " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
        System.out.println();
        
        r1.transform(tensor2);
        System.out.println("r1_afterInvert_andRotate360 " + r1.toString());
        //System.out.println("tensor2 " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
    }
    
}//end of Space3D
