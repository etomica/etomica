package etomica.space.continuum;
import etomica.Space;
import etomica.Constants;
import etomica.Simulation;
    
public class Orientation extends Space.Orientation {
        //The rotation matrix A operates on the components of a vector in the space-fixed frame to yield the
        //components in the body-fixed frame
        private final double[][] A;
    public Orientation() {
        this(3);
    }
    public Orientation(int D) {
        A = new double[D][D];
    }
    private final Vector[] bodyFrame = new Vector[] {new Vector(new double[] {1.0,0.0, 0.0}), new Vector(new double[] {0.0,1.0,0.0}), new Vector(new double[] {0.0,0.0,1.0})};
    private final double[] angle = new double[3];
        private boolean needToUpdateA = true;
        public void E(Space.Orientation o) {E((Orientation)o);}
    public void E(Orientation o) {
        angle[0] = o.angle[0];
        angle[1] = o.angle[1];
        angle[2] = o.angle[2];
        needToUpdateA = true;
    }
        public Space.Vector[] bodyFrame() {return bodyFrame;}
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
    public void randomRotation(double t) {
        rotateBy(0, (2.*Simulation.random.nextDouble()-1.0)*t);
        rotateBy(1, (2.*Simulation.random.nextDouble()-1.0)*t);
        rotateBy(2, (2.*Simulation.random.nextDouble()-1.0)*t);
    }
    private final void updateRotationMatrix() {
        //x-rot
        double theta = angle[0]*(Math.PI / 180);
        double ct = Math.cos(theta);
        double st = Math.sin(theta);
        double Nyx = (double) (A[1][0] * ct + A[2][0] * st);
        double Nyy = (double) (A[1][1] * ct + A[2][1] * st);
        double Nyz = (double) (A[1][2] * ct + A[2][2] * st);
        double Nzx = (double) (A[2][0] * ct - A[1][0] * st);
        double Nzy = (double) (A[2][1] * ct - A[1][1] * st);
        double Nzz = (double) (A[2][2] * ct - A[1][2] * st);
        A[1][0] = Nyx;
        A[1][1] = Nyy;
        A[1][2] = Nyz;
        A[2][0] = Nzx;
        A[2][1] = Nzy;
        A[2][2] = Nzz;
        //y-rot
        theta = angle[1]*(Math.PI / 180);
        ct = Math.cos(theta);
        st = Math.sin(theta);
        double Nxx = (double) (A[0][0] * ct + A[2][0] * st);
        double Nxy = (double) (A[0][1] * ct + A[2][1] * st);
        double Nxz = (double) (A[0][2] * ct + A[2][2] * st);
        Nzx = (double) (A[2][0] * ct - A[0][0] * st);
        Nzy = (double) (A[2][1] * ct - A[0][1] * st);
        Nzz = (double) (A[2][2] * ct - A[0][2] * st);
        A[0][0] = Nxx;
        A[0][1] = Nxy;
        A[0][2] = Nxz;
        A[2][0] = Nzx;
        A[2][1] = Nzy;
        A[2][2] = Nzz;
        //z-rot
        theta = angle[2]*(Math.PI / 180);
        ct = Math.cos(theta);
        st = Math.sin(theta);
        Nyx = (double) (A[1][0] * ct + A[0][0] * st);
        Nyy = (double) (A[1][1] * ct + A[0][1] * st);
        Nyz = (double) (A[1][2] * ct + A[0][2] * st);
        Nxx = (double) (A[0][0] * ct - A[1][0] * st);
        Nxy = (double) (A[0][1] * ct - A[1][1] * st);
        Nxz = (double) (A[0][2] * ct - A[1][2] * st);
        A[1][0] = Nyx;
        A[1][1] = Nyy;
        A[1][2] = Nyz;
        A[0][0] = Nxx;
        A[0][1] = Nxy;
        A[0][2] = Nxz;
        bodyFrame[0].E(A[0]);
        bodyFrame[1].E(A[1]);
        bodyFrame[2].E(A[2]);
        needToUpdateA = false;
        }
    //   public double[][] rotationMatrix() {return A;}
        public void convertToBodyFrame(Vector v) {
/*            if(needToUpdateA) updateRotationMatrix();
        v.x = A[0][0]*v.x + A[0][1]*v.y + A[0][2]*v.z;
        v.y = A[1][0]*v.x + A[1][1]*v.y + A[1][2]*v.z;
        v.z = A[2][0]*v.x + A[2][1]*v.y + A[2][2]*v.z;*/
        }
        public void convertToSpaceFrame(Vector v) {
      /*      if(needToUpdateA) updateRotationMatrix();
        v.x = A[0][0]*v.x + A[1][0]*v.y + A[2][0]*v.z;
        v.y = A[0][1]*v.x + A[1][1]*v.y + A[2][1]*v.z;
        v.z = A[0][2]*v.x + A[1][2]*v.y + A[2][2]*v.z;*/
        }
        public void convertToBodyFrame(Space.Vector v) {convertToBodyFrame((Vector)v);}
        public void convertToSpaceFrame(Space.Vector v) {convertToSpaceFrame((Vector)v);}
}

/*public static class Orientation extends Space.Orientation {
    public final Vector[] eArray = new Vector[] {new Vector(1.0,0.0,0.0)};
    public final Vector L = new Vector(); //angular momentum
    public final Vector e = eArray[0];
    private final Vector[] bodyFrame = new Vector[] {new Vector(1.0,0.0,0.0), 
                                                        new Vector(0.0,1.0,0.0),
                                                        new Vector(0.0,0.0,1.0)};
    public final double[] angle = new double[1];
    public Space.Vector[] direction() {return eArray;}
    public double[] angle() {angle[0] = Math.atan2(e.y,e.x); return angle;}
    public Space.Vector momentum() {return L;}
    public Space.Vector[] bodyFrame() {return bodyFrame;} 
    public void freeFlight(double I) {//motion of a spherical top
    }
}
*/
