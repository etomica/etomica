package etomica.space3d;

import etomica.EtomicaInfo;
import etomica.Space;
import etomica.space.Boundary;

public class Space3D extends Space {

    public final int D() {return D;}
    public final int powerD(int n) {return n*n*n;}
    public final double powerD(double a) {return a*a*a;}
    public int[] makeArrayD(int i) {return new int[] {i, i, i};}
    public double[] makeArrayD(double d) {return new double[] {d, d, d};}
    
    public static final Space3D INSTANCE = new Space3D();
    
    public Space3D() {super(3);}
    
    public double sphereVolume(double r) {return (Math.PI*4.0*r*r*r/3.0);}
    public double sphereArea(double r)  {return (Math.PI*4*r*r);}
    public etomica.space.Vector makeVector() {return new Vector3D();}
    public etomica.space.Orientation makeOrientation() {return new Orientation();}
    public etomica.space.Tensor makeTensor() {return new Tensor3D();}
    public etomica.space.Tensor makeRotationTensor() {return new RotationTensor();}

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Three-dimensional space");
        return info;
    }

    public static final double r2(Vector3D u1, Vector3D u2, Boundary b) {
    	return r2(u1, u2, b, new Vector3D());
    }
    public static final double r2(Vector3D u1, Vector3D u2, Boundary b, Vector3D work) {
    	work.Ev1Mv2(u1, u2);
        b.nearestImage(work);
        return work.squared();
    }  
    
    public static void main (String[] args) {
        Vector3D r1 = new Vector3D(2,2,3);
        System.out.println("r1_before" + r1.toString());
        Tensor3D tensor = new Tensor3D(new double[] {1,2,0,1,1,2,0,0,1});
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
