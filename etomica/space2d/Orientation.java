package etomica.space2d;

import etomica.Constants;
import etomica.Simulation;
import etomica.Space2D;




/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class Orientation extends Orientation {
    //The rotation matrix A operates on the components of a vector in the space-fixed frame to yield the
    //components in the body-fixed frame
    private final double[][] A = new double[Space2D.D][Space2D.D];
    private final Vector[] bodyFrame = new Vector[] {new Vector(1.0,0.0), new Vector(0.0,1.0)};
    private final double[] angle = new double[1];
    private boolean needToUpdateA = true;
    public void E(Orientation o) {E((Orientation)o);}
    public void E(Orientation o) {angle[0] = o.angle[0]; needToUpdateA = true;}
    public Vector[] bodyFrame() {
        if(needToUpdateA) updateRotationMatrix();
        return bodyFrame;
    }
    public double[] angle() {return angle;}
    public final void rotateBy(int i, double dt) {
        if(i == 0) rotateBy(dt);
    }
    public final void rotateBy(double[] dt) {
        rotateBy(dt[0]);
    }
    public final void rotateBy(double dt) {
        angle[0] += dt; 
        if(angle[0] > Constants.TWO_PI) angle[0] -= Constants.TWO_PI;
        else if(angle[0] < 0.0) angle[0] += Constants.TWO_PI;
        needToUpdateA = true;
    }
    public final void randomRotation(double t) {
        rotateBy((2.*Simulation.random.nextDouble()-1.0)*t);
    }
    private final void updateRotationMatrix() {
        A[0][0] = A[1][1] = Math.cos(angle[0]);
        A[0][1] = Math.sin(angle[0]);
        A[1][0] = -A[0][1];
        bodyFrame[0].E(A[0]);
        bodyFrame[1].E(A[1]);
        needToUpdateA = false;
    }
 //   public double[][] rotationMatrix() {return A;}
    public void convertToBodyFrame(Vector v) {
        if(needToUpdateA) updateRotationMatrix();
        double x = A[0][0]*v.x + A[0][1]*v.y;
        v.y = A[1][0]*v.x + A[1][1]*v.y;
        v.x = x;
    }
    public void convertToSpaceFrame(Vector v) {
        if(needToUpdateA) updateRotationMatrix();
        double x = A[0][0]*v.x + A[1][0]*v.y;
        v.y = A[0][1]*v.x + A[1][1]*v.y;
        v.x = x;
    }
    public void convertToBodyFrame(Vector v) {convertToBodyFrame((Vector)v);}
    public void convertToSpaceFrame(Vector v) {convertToSpaceFrame((Vector)v);}
}