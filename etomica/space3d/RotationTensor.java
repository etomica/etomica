package etomica.space3d;

import etomica.Simulation;




/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class RotationTensor extends Tensor3D implements etomica.space.RotationTensor {
    public RotationTensor() {super(); reset();}
    public void reset() {
        xx = 1.0; xy = 0.0; xz = 0.0;
        yx = 0.0; yy = 1.0; yz = 0.0;
        zx = 0.0; zy = 0.0; zz = 1.0;
    }
    /**
     * Sets tensor to for rotation by the given angle about a randomly selected axis.
     */
    public void setAxial(double theta) {
        int n = (int)(Simulation.random.nextDouble()*3);
        setAxial(n, theta);
    }
    /**
     * Sets tensor for rotation about the indicated axis (0=x,1=y,2=z) by 
     * the given angle.
     */
    public void setAxial(int i, double theta) {
        double st = Math.sin(theta);
        double ct = Math.cos(theta);
        switch(i) {
            case 0: xx = 1.; xy = 0.; xz = 0.;
                    yx = 0.; yy = ct; yz = -st;
                    zx = 0.; zy = st; zz = ct;
                    return;
            case 1: xx = ct; xy = 0.; xz = -st;
                    yx = 0.; yy = 1.; yz = 0.;
                    zx = st; zy = 0.; zz = ct;
                    return;
            case 2: xx = ct; xy = -st; xz = 0.;
                    yx = st; yy = ct;  yz = 0.;
                    zx = 0.; zy = 0.;  zz = 1.;
                    return;
            default: throw new IllegalArgumentException("Improper axis specified for Space3D.RotationTensor.setAxial");
        }
    }
    /**
     * Not yet implemented.
     */
    public void setAngles(double[] angles) {
        throw new RuntimeException("Space3D.CoordinateGroup.setAngles() not yet implemented");
    }
    public void invert() {
        double det = xx*yy*zz - xx*yz*zy - yx*xy*zz + yx*xz*zy + zx*xy*yz - zx*xz*yy ;
        double xx1 = (yy*zz - yz*zy)/det;
        double xy1 = (-xy*zz + xz*zy)/det;
        double xz1 = (xy*yz - xz*yy)/det;
        
        double yx1 = (-yx*zz + yz*zx)/det;
        double yy1 = (xx*zz - xz*zx)/det;
        double yz1 = (-xx*yz + xz*yx)/det;
        
        double zx1 = (yx*zy - yy*zx)/det;
        double zy1 = (-xx*zy + xy*zx)/det;
        double zz1 = (xx*yy - xy*yx)/det;
        
        this.xx = xx1; this.xy = xy1; this.xz = xz1;
        this.yx = yx1; this.yy = yy1; this.yz = yz1;
        this.zx = zx1; this.zy = zy1; this.zz = zz1;
        
    }
    /**
     * Method to test rotation tensor.
     */
    public static void main (String[] args) {
        
        Vector3D r1 = new Vector3D(2,2,3);
        System.out.println("r1_before" + r1.toString());
        Tensor3D tensor = new Tensor3D(new double[] {1,2,0,1,1,2,0,0,1});
        RotationTensor tensor2 = new RotationTensor();
        tensor2.E(tensor);
        System.out.println("tensor2_before " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
        System.out.println();
    
        r1.transform(tensor2);
        System.out.println("r1_transform(tensor2)" + r1.toString());
        tensor2.invert();
        System.out.println("tensor2_invert " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
        //tensor2.setAxial(1, 2*Math.PI);
        //System.out.println("tensor2_rotate_360 " + tensor2.xx + "  " +tensor2.xy +"  "+tensor2.xz +"  "+tensor2.yx +"  "+tensor2.yy +"  "+tensor2.yz +"  "+tensor2.zx +"  "+tensor2.zy +"  "+tensor2.zz); 
        //System.out.println();
    
        //r1.transform(tensor2);
        //System.out.println("r1_afterInvert_andRotate360 " + r1.toString());
    }//end of main 
    
}
