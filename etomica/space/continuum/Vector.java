package etomica.space.continuum;
import etomica.Space;
import etomica.Space3D;
import etomica.Space2D;
import etomica.Simulation;

public final class Vector extends Space.Vector {    

    public final double x[];
    private final int D;

    public static final Vector WORK = new Vector();

    public Vector () {this(3);}
    
    public Vector(int d) {
        x = new double[d];
        D = d;
    }

    public Vector (double[] a) {
        this(a.length);
        for(int i=0; i<a.length; i++) x[i] = a[i];
    }

    public Vector (Space.Vector u) {
        this(u.toArray());
    }

    public double[] toArray() {return x;}

    public void sphericalCoordinates(double[] result) {
        System.out.println("Vector: Shouldn't be here");
    }

    public int length() {return D;}

    public int D() {return D;}

    public double component(int i) {return x[i];}

    public void setComponent(int i, double d) {x[i]=d;}

    public void E(Space.Vector u) {
        double[] ux;
        try {
          ux = ((Vector)u).x;
        } catch(ClassCastException e) {ux = u.toArray();}
        switch(D) {
            case 3: x[2] = ux[2];
            case 2: x[1] = ux[1];
            case 1: x[0] = ux[0];
        }
    }
       
    public void E(double[] u) {//should check length of array for exception
        switch(D) {
            case 3: x[2] = u[2];
            case 2: x[1] = u[1];
            case 1: x[0] = u[0];
        }
    }  

    public void E(double a) {
        switch(D) {
            case 3: x[2] = a;
            case 2: x[1] = a;
            case 1: x[0] = a;
        }
    }

    public void Ea1Tv1(double a1, Space.Vector u) {
        double[] ux = ((Vector)u).x;
        switch(D) {
            case 3: x[2] = a1*ux[2];
            case 2: x[1] = a1*ux[1];
            case 1: x[0] = a1*ux[0];
        }
    }

    public void PEa1Tv1(double a1, Space.Vector u) {
        double[] ux = ((Vector)u).x;
        switch(D) {
            case 3: x[2] += a1*ux[2];
            case 2: x[1] += a1*ux[1];
            case 1: x[0] += a1*ux[0];
        }
    }

    public void PE(Space.Vector u) {
        double[] ux = ((Vector)u).x;
        switch(D) {
            case 3: x[2] += ux[2];
            case 2: x[1] += ux[1];
            case 1: x[0] += ux[0];
        }
    }

    public void PE(double a) {
  /*          x[2] += a;
            x[1] += a;
            x[0] += a;
*/        switch(D) {
            case 3: x[2] += a;
            case 2: x[1] += a;
            case 1: x[0] += a;
        }
    }

    public void ME(Space.Vector u) {
        double[] ux = ((Vector)u).x;
 /*           x[2] -= ux[2];
            x[1] -= ux[1];
            x[0] -= ux[0];
*/        switch(D) {
            case 3: x[2] -= ux[2];
            case 2: x[1] -= ux[1];
            case 1: x[0] -= ux[0];
        }
    }

    public void PE(int i, double a) {x[i]+=a;}

    public void TE(double a) {
        switch(D) {
            case 3: x[2] *= a;
            case 2: x[1] *= a;
            case 1: x[0] *= a;
        }
    }

    public void TE(Space.Vector u) {
        double[] ux = ((Vector)u).x;
        switch(D) {
            case 3: x[2] *= ux[2];
            case 2: x[1] *= ux[1];
            case 1: x[0] *= ux[0];
        }
    }

    public void TE(int i, double a) {x[i]*=a;}

    //get rid of this
    public void DE(double a) {for(int i=D()-1;i>=0;i--)x[i]/=a;}

    public void DE(Space.Vector u) {
        double[] ux = ((Vector)u).x;
        switch(D) {
            case 3: x[2] /= ux[2];
            case 2: x[1] /= ux[1];
            case 1: x[0] /= ux[0];
        }
    }
    
    public void Ev1Pv2(Space.Vector u1, Space.Vector u2) {
        double[] x1 = ((Vector)u1).x;
        double[] x2 = ((Vector)u2).x;
        switch(D) {
            case 3: x[2] = x1[2] + x2[2];
            case 2: x[1] = x1[1] + x2[1];
            case 1: x[0] = x1[0] + x2[0];
        }
    }

    public void Ev1Mv2(Space.Vector u1, Space.Vector u2) {
        double[] x1 = ((Vector)u1).x;
        double[] x2 = ((Vector)u2).x;
/*            x[2] = x1[2] - x2[2];
            x[1] = x1[1] - x2[1];
            x[0] = x1[0] - x2[0];
*/        switch(D) {
            case 3: x[2] = x1[2] - x2[2];
            case 2: x[1] = x1[1] - x2[1];
            case 1: x[0] = x1[0] - x2[0];
        }
    }

    public double min() {
        double m = Double.MAX_VALUE;
        switch(D) {
            case 3: m = x[2];
            case 2: m = (m < x[1]) ? m : x[1];
            case 1: m = (m < x[0]) ? m : x[0];
        }
        return m;
    }

    public double max() {
        double m = -Double.MAX_VALUE;
        switch(D) {
            case 3: m = x[2];
            case 2: m = (m > x[1]) ? m : x[1];
            case 1: m = (m > x[0]) ? m : x[0];
        }
        return m;
    }
    
    public double squared() {
        switch(D) {
            case 3: return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
            case 2: return x[0]*x[0] + x[1]*x[1];
            case 1: return x[0]*x[0];
            default: throw new RuntimeException();
        }
    }
    
    public void mod(double a) {
        switch(D) {
            case 3:
                while(x[2] > a)   x[2] -= a;
                while(x[2] < 0.0) x[2] += a;
            case 2:
                while(x[1] > a)   x[1] -= a;
                while(x[1] < 0.0) x[1] += a;
            case 1:
                while(x[0] > a)   x[0] -= a;
                while(x[0] < 0.0) x[0] += a;
        }
    }

    public void mod(Space.Vector u) {
        double[] ux = ((Vector)u).x;
        double[] x = this.x;
        double uu;
        switch(D) {
            case 3:
                uu = ux[2];
                while(x[2] > uu) x[2] -= uu;
                while(x[2] < 0.) x[2] += uu;
            case 2:
                uu = ux[1];
                while(x[1] > uu) x[1] -= uu;
                while(x[1] < 0.) x[1] += uu;
            case 1:
                uu = ux[0];
                while(x[0] > uu) x[0] -= uu;
                while(x[0] < 0.) x[0] += uu;
        }
    }
    
    public double dot(Space.Vector u) {
        double[] ux = ((Vector)u).x;
        switch(D) {
            case 3: return x[0]*ux[0] + x[1]*ux[1] + x[2]*ux[2];
            case 2: return x[0]*ux[0] + x[1]*ux[1];
            case 1: return x[0]*ux[0];
            default: throw new RuntimeException();
        }
    }
    
    public double dot(double ux[]){
        switch(D) {
            case 3: return x[0]*ux[0] + x[1]*ux[1] + x[2]*ux[2];
            case 2: return x[0]*ux[0] + x[1]*ux[1];
            case 1: return x[0]*ux[0];
            default: throw new RuntimeException();
        }
    }
    public Space.Vector P(Space.Vector u) {
        double[] Wx = WORK.x;
        double[] ux = ((Vector)u).x;
        switch(D) {
            case 3: Wx[2] = x[2] + ux[2];
            case 2: Wx[1] = x[1] + ux[1];
            case 1: Wx[0] = x[0] + ux[0];
        }
        return WORK;
    }

    public Space.Vector M(Space.Vector u) {
        double[] Wx = WORK.x;
        double[] ux = ((Vector)u).x;
        switch(D) {
            case 3: Wx[2] = x[2] - ux[2];
            case 2: Wx[1] = x[1] - ux[1];
            case 1: Wx[0] = x[0] - ux[0];
        }
        return WORK;
    }
        

    public Space.Vector T(Space.Vector u) {
        double[] Wx = WORK.x;
        double[] ux = ((Vector)u).x;
        switch(D) {
            case 3: Wx[2] = x[2] * ux[2];
            case 2: Wx[1] = x[1] * ux[1];
            case 1: Wx[0] = x[0] * ux[0];
        }
        return WORK;
    }

    public Space.Vector D(Space.Vector u) {
        double[] Wx = WORK.x;
        double[] ux = ((Vector)u).x;
        switch(D) {
            case 3: Wx[2] = x[2] / ux[2];
            case 2: Wx[1] = x[1] / ux[1];
            case 1: Wx[0] = x[0] / ux[0];
        }
        return WORK;
    }
        

    public Space.Vector abs() {
        double[] Wx = WORK.x;
        switch(D) {
            case 3: Wx[2] = (x[2] > 0) ? x[2] : -x[2];
            case 2: Wx[1] = (x[1] > 0) ? x[1] : -x[1];
            case 1: Wx[0] = (x[0] > 0) ? x[0] : -x[0];
        }
        return WORK;
    }

    /**
     * Sets vector to a random point in a unit cube, centered on the origin.
     */
    public void setRandomCube(){
        switch(D) {
            case 3: 
                x[0] = Simulation.random.nextDouble()-0.5;
                x[1] = Simulation.random.nextDouble()-0.5;
                x[2] = Simulation.random.nextDouble()-0.5;
                break;
            case 2:
                x[0] = Simulation.random.nextDouble()-0.5;
                x[1] = Simulation.random.nextDouble()-0.5;
                break;
            case 1:
                x[0] = Simulation.random.nextDouble()-0.5;
                break;
        }
    }

    /**
     * Replaces this vector with its cross-product with the given 3D vector, with result projected
     * onto the 2D plane.  This vector becomes the result of (this vector) x u.
     */

    public void XE(Space3D.Vector u) {
        throw new RuntimeException("Method Vector.XE not implemented");
    }

    public Space3D.Vector cross(Space3D.Vector u) {//not thread safe
        throw new RuntimeException("Method Vector.cross not implemented");
    }
        /*
        Space3D.Vector.WORK.x = y*u.z;
        Space3D.Vector.WORK.y = -x*u.z;
        Space3D.Vector.WORK.z = x*u.y - y*u.x;
        return Space3D.Vector.WORK;
        */
    public Space3D.Vector cross(Space2D.Vector u) {//not thread safe
        throw new RuntimeException("Method Vector.cross not implemented");
    }

    public void normalize() {
        this.TE(1.0/Math.sqrt(this.squared()));
    }

    public void transform(Space.Tensor A) {transform((Space.Tensor)A);}

    public void transform(Space.Boundary b, Space.Vector r0, Space.Tensor A) {
        System.out.println("Vector: shouldn't be here 0");
    }

    public void randomStep(double d) {System.out.println("Vector: shouldn't be here 1");} //uniformly distributed random step in x and y, within +/- d

    public void setRandom(double d) {System.out.println("Vector: shouldn't be here 2");; }

    public void setRandom(double dx, double dy) {System.out.println("Vector: shouldn't be here 3");;}

    public void setRandom(Vector u) {System.out.println("Vector: shouldn't be here 4");;}

    /**
     * Sets vector to a random point within a sphere of unit diameter, centered on the origin.
     */
    public void setRandomSphere() {
        switch(D) {
            case 1:
                x[0] = Simulation.random.nextDouble() - 0.5;
                break;
            case 2:
                double z1 = 0.0;
                double z2 = 0.0;
                double rsq = Double.MAX_VALUE;
                while(rsq > 0.25) {
                    z1 = Simulation.random.nextDouble() - 0.5;
                    z2 = Simulation.random.nextDouble() - 0.5;
                    rsq = z1*z1 + z2*z2;
                }
                x[0] = z1;;
                x[1] = z2;;
                break;
            case 3:
                double y1 = 0.0;
                double y2 = 0.0;
                double y3 = 0.0;
                double r2 = Double.MAX_VALUE;
                while(r2 > 0.25) {
                    y1 = Simulation.random.nextDouble() - 0.5;
                    y2 = Simulation.random.nextDouble() - 0.5;
                    y3 = Simulation.random.nextDouble() - 0.5;
                    r2 = y1*y1 + y2*y2 + y3*y3;
                }
                x[0] = y1;
                x[1] = y2;;
                x[2] = y3;;
        }
        //compact but less efficient formulation
            //do this.setRandomCube();
            //while(this.squared() > 0.25);
    }//end of setRandomSphere

    public void randomRotate(double thetaStep){
        System.out.println("Vector: shouldn't be here 6");
        }

        

    public double[] getArray(){return x;}
}

