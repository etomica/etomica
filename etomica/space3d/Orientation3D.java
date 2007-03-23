package etomica.space3d;

import etomica.space.IVector;
import etomica.util.Constants;
import etomica.util.IRandom;

public class Orientation3D extends etomica.space.Orientation {
    public Orientation3D() {
        throw new RuntimeException("Space3D.Orientation should be checked for correctness before using");
    }

    private static final long serialVersionUID = 1L;
    //The rotation matrix A operates on the components of a vector in the space-fixed frame to yield the
    //components in the body-fixed frame
    private final double[][] A = new double[3][3];
    private final Vector3D[] bodyFrame = new Vector3D[] {new Vector3D(1.0,0.0, 0.0), new Vector3D(0.0,1.0,0.0), new Vector3D(0.0,0.0,1.0)};
    private final double[] angle = new double[3];
    private final Vector3D orientVector = new Vector3D(1.0,0.0,0.0);
    private boolean needToUpdateA = true;
    private transient double x1,y1,z1;//temp variable
    private transient Vector3D v1 = new Vector3D();
    public void E(etomica.space.Orientation o) {E((Orientation3D)o);}
    public IVector getOrientation(){return orientVector;}
    public void setOrientation(Vector3D vect){orientVector.E(vect);}
    public void E(Orientation3D o) {
      angle[0] = o.angle[0];
      angle[1] = o.angle[1];
      angle[2] = o.angle[2];
      needToUpdateA = true;
      orientVector.E(o.getOrientation());
    }
    public IVector[] bodyFrame() {return bodyFrame;}
    public double[] angle() {return angle;}
    public final void rotateBy(double dt[]) {
        rotateBy(0, dt[0]);
        rotateBy(1, dt[1]);
        rotateBy(2, dt[2]);
    }
    public final void rotateBy(double dt) {
      rotateBy(0, dt);
      rotateBy(1, dt);
      rotateBy(2, dt);
    }
    public void rotateBy(int i, double dt) {
        angle[i] += dt;
        if(angle[i] > Constants.TWO_PI) angle[i] -= Constants.TWO_PI;
        else if(angle[i] < 0.0) angle[i] += Constants.TWO_PI;
        needToUpdateA = true;
    }
      /**
    * Choose one of the three spaced fixed axes at random and rotate 
    * it by t (in radian).
    */
    public void randomRotation(IRandom random, double stepsize) {
        double t = (2.0*random.nextDouble()-1.0)*stepsize;
        int i = (int)(3.0 *(random.nextDouble()));
        double ct,st;
        
        ct= Math.cos(t);st =Math.sin(t);
        
        //rotateBy(i,t);
        
        switch(i){
         case 0:
                x1=orientVector.x(0);
                y1=ct*orientVector.x(1)+ st*orientVector.x(2);      
                z1=ct*orientVector.x(2)-st*orientVector.x(1);
                break;
         case 1:
                x1=ct*orientVector.x(0)-st*orientVector.x(2);
                y1=orientVector.x(1);
                z1=ct*orientVector.x(2)+st*orientVector.x(0);
                break;
         case 2:
                x1=ct*orientVector.x(0)+st*orientVector.x(1);
                y1=ct*orientVector.x(1)-st*orientVector.x(0);
                z1=orientVector.x(2);
                break;
        }
            
        orientVector.E(x1,y1,z1);
       // orientVector.normalize(); // Not needed
        
            
            
        
    }
    private final void updateRotationMatrix() {
        double ct0,ct1,ct2;
        double st0,st1,st2;
                    
        ct0=Math.cos(angle[0]);ct1=Math.cos(angle[1]);ct2=Math.cos(angle[2]);
        st0=Math.sin(angle[0]);st1=Math.sin(angle[1]);st2=Math.sin(angle[2]);
                 
        //LEFT HAND SYSTEM  
        A[0][0]= ct1*ct2;
        A[0][1]= -ct1*st2;
        A[0][2]= st1;
        A[1][0]= st0*st1*ct2+ct0*st2;
        A[1][1]= -st0*st1*st2+ct0*ct2;
        A[1][2]= -st0*ct1;
        A[2][0]= -ct0*st1*ct2+st0*st2;
        A[2][1]= ct0*st1*st2+st0*ct2;
        A[2][2]= ct0*ct1;
                  
        bodyFrame[0].E(A[0]);
        bodyFrame[1].E(A[1]);
        bodyFrame[2].E(A[2]);
        needToUpdateA = false;
     }
  //   public double[][] rotationMatrix() {return A;}
      public void convertToBodyFrame(IVector[] u) {
        if(needToUpdateA) updateRotationMatrix();
        Vector3D[] v = (Vector3D[])u;
        for(int i=0; i<v.length; i++) {
            v1.x = A[0][0]*v[i].x + A[0][1]*v[i].y + A[0][2]*v[i].z;
            v1.y = A[1][0]*v[i].x + A[1][1]*v[i].y + A[1][2]*v[i].z;
            v1.z = A[2][0]*v[i].x + A[2][1]*v[i].y + A[2][2]*v[i].z;
            v[i].E(v1);
        }
     }
     //V_space = A_transpose*V_body
     public void convertToSpaceFrame(IVector[] u) {
        if(needToUpdateA) updateRotationMatrix();
        Vector3D[] v = (Vector3D[])u;
        for(int i=0; i<v.length; i++) {
            v1.x = A[0][0]*v[i].x + A[1][0]*v[i].y + A[2][0]*v[i].z;
            v1.y = A[0][1]*v[i].x + A[1][1]*v[i].y + A[2][1]*v[i].z;
            v1.z = A[0][2]*v[i].x + A[1][2]*v[i].y + A[2][2]*v[i].z;
            v[i].E(v1);
        }
     }
    public void convertToBodyFrame(IVector v) {convertToBodyFrame(v);}
    public void convertToSpaceFrame(IVector v) {convertToSpaceFrame(v);}
}
