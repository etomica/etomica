package etomica.space.continuum;
import etomica.Space;
import etomica.Space3D;
import etomica.Space2D;
import etomica.Simulation;

//not completed/debugged

public final class Vector3D extends Space.Vector {    

    public final double X[];
    public double x, y, z;
    private final int D;

    public static final Vector WORK = new Vector();

    public Vector3D () {this(3);}
    
    public Vector3D(int d) {
        X = new double[d];
        D = d;
    }

    public Vector3D (double[] a) {
        this(a.length);
        x = a[0];
        y = a[1];
        z = a[2];
    }

    public Vector3D(Space.Vector u) {
        this(u.toArray());
    }
    
    private void copy() {
        X[0] = x;
        X[1] = y;
        X[2] = z;
    }
    private void uncopy() {
        x = X[0];
        y = X[1];
        z = X[2];
    }

    public double[] toArray() {copy(); return X;}

    public void sphericalCoordinates(double[] result) {
        System.out.println("Vector: Shouldn't be here");
    }

    public int length() {return D;}//x.length;}

    public int D() {return D;}//x.length;}

    public double component(int i) {copy(); return X[i];}

    public void setComponent(int i, double d) {copy(); X[i]=d; uncopy();}

    public void E(Space.Vector u1) {
        Vector3D u;
        try{ u = (Vector3D)u1;}
        catch(ClassCastException e) {this.E((Space3D.Vector)u1); return;}
        x = u.x;
        y = u.y;
        z = u.z;
    }
    public void E(Space3D.Vector u) {
        x = u.x;
        y = u.y;
        z = u.z;
    }
       
    public void E(double[] u) {
        x=u[0];
        y=u[1];
        z=u[2];
    }  //should check length of array for exception

    public void E(double a) {x=y=z=a;}//for(int i=D()-1;i>=0;i--)x[i]=a;}

    public void Ea1Tv1(double a1, Space.Vector u) {
        Vector3D u1 = (Vector3D)u;
        x = a1*u1.x;
        y = a1*u1.y;
        z = a1*u1.z;
    }

    public void PEa1Tv1(double a1, Space.Vector u) {
        Vector3D u1 = (Vector3D)u;
        x += a1*u1.x;
        y += a1*u1.y;
        z += a1*u1.z;
    }

    public void PE(Space.Vector u1) {
        Vector3D u = (Vector3D)u1;
        x += u.x;
        y += u.y;
        z += u.z;
    }

    public void PE(double a) {
        x += a;
        y += a;
        z += a;
    }

    public void ME(Space.Vector u1) {
        Vector3D u = (Vector3D)u1;
        x -= u.x;
        y -= u.y;
        z -= u.z;
    }

    public void PE(int i, double a) {}//x[i]+=a;}

    public void TE(double a) {}//for(int i=D()-1;i>=0;i--)x[i]*=a;}

    public void TE(Space.Vector u) {}//for(int i=D()-1;i>=0;i--)x[i]*=u.component(i);}

        

    public void TE(int i, double a) {}//x[i]*=a;}

    //get rid of this
    public void DE(double a) {}//for(int i=D()-1;i>=0;i--)x[i]/=a;}

    public void DE(Space.Vector u) {}//for(int i=D()-1;i>=0;i--)x[i]/=u.component(i);}
    
    public void Ev1Pv2(Space.Vector u1, Space.Vector u2) {
        //for(int i=D()-1; i>=0; i--) x[i] = u1.component(i) + u2.component(i);
    }

    public void Ev1Mv2(Space.Vector u11, Space.Vector u22) {
        Vector3D u1 = (Vector3D)u11;
        Vector3D u2 = (Vector3D)u22;
        x = u1.x - u2.x;
        y = u1.y - u2.y;
        z = u1.z - u2.z;
    }

    public double min() {return 0.0;}//double m=Double.MAX_VALUE;for(int i=D()-1;i>=0;i--)if(x[i]<m) m=x[i]; return m;}

    public double max() {return 0.0;}//double m=-Double.MAX_VALUE;for(int i=D()-1;i>=0;i--)if(x[i]>m) m=x[i]; return m;}

    public double squared() {
            return x*x + y*y + z*z;
    }
    
    public void mod(double a) {
            while(z > a)   z -= a;
            while(z < 0.0) z += a;
            while(y > a)   y -= a;
            while(y < 0.0) y += a;
            while(x > a)   x -= a;
            while(x < 0.0) x += a;
    }

    public void mod(Space.Vector u1) {
        /*double uu;
        Vector3D u = (Vector3D)u1;
        uu = u.z;
        while(z > uu) z -= uu;
        while(z < 0.) z += uu;
        uu = u.y;
        while(y > uu) y -= uu;
        while(y < 0.) y += uu;
        uu = u.x;
        while(x > uu) x -= uu;
        while(x < 0.) x += uu;*/
        Vector3D u = (Vector3D)u1;
        while(z > u.z) z -= u.z;
        while(z < 0.0) z += u.z;
        while(y > u.x) y -= u.y;
        while(y < 0.0) y += u.y;
        while(x > u.x) x -= u.x;
        while(x < 0.0) x += u.x;
    }
    
//    public void mod(Space.Vector u1, Space.Vector v1) {
//        Vector3D u = (Vector3D)u1;
//        Vector3D v = (Vector3D)v1;
    public void mod(Vector3D u, Vector3D v) {
        while(z > +u.z) z -= v.z;
        while(z < -u.z) z += v.z;
        while(y > +u.y) y -= v.y;
        while(y < -u.y) y += v.y;
        while(x > +u.x) x -= v.x;
        while(x < -u.x) x += v.x;
    }

    public double dot(Space.Vector u1) {
        Vector3D u = (Vector3D)u1;
        return x*u.x + y*u.y + z*u.z;
    }
    public double dot(double a[]){return x*a[0]+y*a[1]+z*a[2];}
    public Space.Vector P(Space.Vector u) {return WORK;}//NEED TO CHANGE

    public Space.Vector M(Space.Vector u) { return WORK;}//NEED TO CHANGE

    public Space.Vector T(Space.Vector u) { return WORK;}//NEED TO CHANGE

    public Space.Vector D(Space.Vector u) {return WORK;} //NEED TO CHANGE

    public Space.Vector abs() { System.out.println("Vector.abs shouldn't be here");return WORK;} // NEED TO CHANGE

    public void setRandomCube(){
        x = Simulation.random.nextDouble()-0.5;
        y = Simulation.random.nextDouble()-0.5;
        z = Simulation.random.nextDouble()-0.5;
    }

    /**
     * Replaces this vector with its cross-product with the given 3D vector, with result projected
     * onto the 2D plane.  This vector becomes the result of (this vector) x u.
     */

     public void XE(Space3D.Vector u) {
     }

    //NEED TO CHANGE

       public Space3D.Vector cross(Space3D.Vector u) {//not thread safe
        return null;
       }
       public Space3D.Vector cross(Space2D.Vector u) {//not thread safe
        return null;
       }
        /*

        Space3D.Vector.WORK.x = y*u.z;

        Space3D.Vector.WORK.y = -x*u.z;

        Space3D.Vector.WORK.z = x*u.y - y*u.x;

        */

//           return Space3D.Vector.WORK;

//       }

    //NEED TO CHANGE

    public void normalize() {
    }

    public void transform(Space.Tensor A) {transform((Space.Tensor)A);}

    public void transform(Space.Boundary b, Space.Vector r0, Space.Tensor A) {
        System.out.println("Vector: shouldn't be here 0");
    }

    public void randomStep(double d) {System.out.println("Vector: shouldn't be here 1");} //uniformly distributed random step in x and y, within +/- d

    public void setRandom(double d) {System.out.println("Vector: shouldn't be here 2");; }

    public void setRandom(double dx, double dy) {System.out.println("Vector: shouldn't be here 3");;}

    public void setRandom(Vector u) {System.out.println("Vector: shouldn't be here 4");;}

    public void setRandomSphere() {///check this method
            double z1 = 0.0;
            double z2 = 0.0;
            double z3 = 0.0;
            double rsq = Double.MAX_VALUE;
            while(rsq > 1.0) {
                
                z1 = 1.0 - 2.0*Simulation.random.nextDouble();
                z2 = 1.0 - 2.0*Simulation.random.nextDouble();
                z3 = 1.0 - 2.0*Simulation.random.nextDouble();
        
                rsq = z1*z1+z2*z2+z3*z3;
            }
            double r = Math.sqrt(rsq);
            x = z1/r;
            y = z2/r;
            z = z3/r;
        }

    public void randomRotate(double thetaStep){
        System.out.println("Vector: shouldn't be here 6");
        }

        

    public double[] getArray(){copy();return X;}
}

