package etomica.space.continuum;
import etomica.Space;
import etomica.Space3D;
import etomica.Space2D;
import etomica.Simulation;

public final class VectorSave extends Space.Vector {    

    public final double x[];
    private final int D;

    public static final Vector WORK = new Vector();

    public VectorSave () {this(3);}
    
    public VectorSave(int d) {
        x = new double[d];
        D = d;
    }

    public VectorSave (double[] a) {
        this(a.length);
        for(int i=0; i<a.length; i++) x[i] = a[i];
    }

    public VectorSave (Vector u) {
        this(u.toArray());
    }

    public double[] toArray() {return x;}

    public void sphericalCoordinates(double[] result) {
        System.out.println("Vector: Shouldn't be here");
    }

    public int length() {return D;}//x.length;}

    public int D() {return D;}//x.length;}

    public double component(int i) {return x[i];}

    public void setComponent(int i, double d) {x[i]=d;}

    public void E(Space.Vector u) {
        x[0] = u.component(0);
        x[1] = u.component(1);
        x[2] = u.component(2);
        //   for(int i=D()-1;i>=0;i--)x[i]=u.x[i];}
    }
       
    public void E(double[] u) {
        x[0]=u[0];
        x[1]=u[1];
        x[2]=u[2];
    }  //should check length of array for exception

    public void E(double a) {x[0]=x[1]=x[2]=a;}//for(int i=D()-1;i>=0;i--)x[i]=a;}

/*        public void E(Space.Vector u) {
            if (u instanceof Space3D.Vector){
            x=u.toArray();
        }
    else 
            E((Vector)u);
    }*/

    public void Ea1Tv1(double a1, Space.Vector u) {
        double[] ux = ((VectorSave)u).x;
        x[0] = a1*ux[0];
        x[1] = a1*ux[1];
        x[2] = a1*ux[2];
        // for(int i=D()-1;i>=0;i--)x[i]=a1*u1.x[i];
        }

    public void PEa1Tv1(double a1, Space.Vector u) {
        double[] ux = ((VectorSave)u).x;
        x[0] += a1*ux[0];
        x[1] += a1*ux[1];
        x[2] += a1*ux[2];
            
    //    for(int i=D()-1;i>=0;i--)x[i]+=a1*u1.x[i];
        }

    public void PE(Space.Vector u) {
        double[] ux = ((VectorSave)u).x;
        x[0] += ux[0];
        x[1] += ux[1];
        x[2] += ux[2];
        //for(int i=D()-1;i>=0;i--)x[i]+=u.x[i];
        }

    public void PE(double a) {
        x[0] += a;
        x[1] += a;
        x[2] += a;
   //     for(int i=D()-1;i>=0;i--)x[i]+=a;
    }

    public void ME(Space.Vector u) {
        x[0] -= u.component(0);
        x[1] -= u.component(1);
        x[2] -= u.component(2);
    //    for(int i=D()-1;i>=0;i--)x[i]-=u.component(i);
    }

    public void PE(int i, double a) {x[i]+=a;}

    public void TE(double a) {for(int i=D()-1;i>=0;i--)x[i]*=a;}

    public void TE(Space.Vector u) {for(int i=D()-1;i>=0;i--)x[i]*=u.component(i);}

        

    public void TE(int i, double a) {x[i]*=a;}

    //get rid of this
    public void DE(double a) {for(int i=D()-1;i>=0;i--)x[i]/=a;}

    public void DE(Space.Vector u) {for(int i=D()-1;i>=0;i--)x[i]/=u.component(i);}
    
    public void Ev1Pv2(Space.Vector u1, Space.Vector u2) {
        for(int i=D()-1; i>=0; i--) x[i] = u1.component(i) + u2.component(i);
    }

    public void Ev1Mv2(Space.Vector u1, Space.Vector u2) {
        double[] x1 = ((VectorSave)u1).x;
        double[] x2 = ((VectorSave)u2).x;
        x[0] = x1[0] - x2[0];
        x[1] = x1[1] - x2[1];
        x[2] = x1[2] - x2[2];
     /*   x[0] = u1.component(0) - u2.component(0);
        x[1] = u1.component(1) - u2.component(1);
        x[2] = u1.component(2) - u2.component(2);
  *///      for(int i=D()-1; i>=0; i--) x[i] = u1.component(i) - u2.component(i);
    }

    public double min() {double m=Double.MAX_VALUE;for(int i=D()-1;i>=0;i--)if(x[i]<m) m=x[i]; return m;}

    public double max() {double m=-Double.MAX_VALUE;for(int i=D()-1;i>=0;i--)if(x[i]>m) m=x[i]; return m;}

    public double squared() {
            return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
   /*     switch(D()) {
            case 1: return x[0]*x[0];
            case 2: return x[0]*x[0] + x[1]*x[1];
            case 3: return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
            default: throw new RuntimeException();
        //   double sum=0.0; for(int i=D()-1;i>=0;i--)sum+=x[i]*x[i]; return sum;
        }
 */   }
    
    public void mod(double a) {
            while(x[2] > a)   x[2] -= a;
            while(x[2] < 0.0) x[2] += a;
            while(x[1] > a)   x[1] -= a;
            while(x[1] < 0.0) x[1] += a;
            while(x[0] > a)   x[0] -= a;
            while(x[0] < 0.0) x[0] += a;
 /*       for(int i=D()-1; i>=0; i--) {
            while(x[i] > a)   x[i] -= a;
            while(x[i] < 0.0) x[i] += a;
        }*/
    }

    public void mod(Space.Vector u) {
        double[] ux = ((VectorSave)u).x;
            double uu = ux[2];
            while(x[2] > uu) x[2] -= uu;
            while(x[2] < 0.) x[2] += uu;
            uu = ux[1];
            while(x[1] > uu) x[1] -= uu;
            while(x[1] < 0.) x[1] += uu;
            uu = ux[0];
            while(x[0] > uu) x[0] -= uu;
            while(x[0] < 0.) x[0] += uu;
/*        switch(D()) {
            case 3: uu = u.component(2);
                    while(x[2] > uu) x[2] -= uu;
                    while(x[2] < 0.) x[2] += uu;
            case 2: uu = u.component(1);
                    while(x[1] > uu) x[1] -= uu;
                    while(x[1] < 0.) x[1] += uu;
            case 1: uu = u.component(0);
                    while(x[0] > uu) x[0] -= uu;
                    while(x[0] < 0.) x[0] += uu;
                    break;
            default: throw new RuntimeException();
        }*/
   //     for(int i=D()-1; i>=0; i--) {
   //         while(x[i] > u.component(i)) x[i] -= u.component(i);
   //         while(x[i] < 0.0)            x[i] += u.component(i);
   //     }
    }

    public double dot(Space.Vector u) {
        return x[0]*u.component(0) + x[1]*u.component(1) + x[2]*u.component(2);
//            double sum=0.0;for(int i=D()-1;i>=0;i--)sum+=x[i]*u.x[i]; return sum;}
    }
    public double dot(double a[]){return x[0]*a[0]+x[1]*a[1]+x[2]*a[2];}
    public Space.Vector P(Space.Vector u) {return WORK;}//NEED TO CHANGE

    public Space.Vector M(Space.Vector u) { return WORK;}//NEED TO CHANGE

    public Space.Vector T(Space.Vector u) { return WORK;}//NEED TO CHANGE

    public Space.Vector D(Space.Vector u) {return WORK;} //NEED TO CHANGE

    public Space.Vector abs() { System.out.println("shouldn't be here");return WORK;} // NEED TO CHANGE

    public void setRandomCube(){
        x[0] = Simulation.random.nextDouble()-0.5;
        x[1] = Simulation.random.nextDouble()-0.5;
        x[2] = Simulation.random.nextDouble()-0.5;
  //      for(int i=D()-1;i>=0;i--) x[i]=Simulation.random.nextDouble()-0.5;
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
            x[0] = z1/r;
            x[1] = z2/r;
            x[2] = z3/r;
        }

    public void randomRotate(double thetaStep){
        System.out.println("Vector: shouldn't be here 6");
        }

        

    public double[] getArray(){return x;}
}

