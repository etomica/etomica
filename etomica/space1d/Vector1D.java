package etomica.space1d;

import etomica.Simulation;
import etomica.space.Boundary;
import etomica.space.Vector;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public final class Vector1D extends etomica.space.Vector {  //declared final for efficient method calls
        public static final Vector1D ZERO = new Vector1D(0.0);
        public static final Vector1D WORK = new Vector1D();
        double x;
        public Vector1D () {x = 0.0;}
        public Vector1D (double a1) {x = a1;}
        public Vector1D (double[] a) {x = a[0];}//should check length of a for exception
        public Vector1D (Vector1D u) {this.E(u);}
        public boolean equals(etomica.space.Vector v) {return equals((Vector1D)v);}
        public boolean equals(Vector1D v) {return (x == v.x);}
        public boolean isZero() {return x==0;}
        public String toString() {return "("+x+")";}
        public int D() {return 1;}
        public double x(int i) {return x;}
        public void setX(int i, double d) {x=d;}
        public void assignTo(double[] array) {array[0] = x;}
        public double[] toArray() {return new double[] {x};}
        public void sphericalCoordinates(double[] result) {result[0] = x;}
        public void E(Vector1D u) {x = u.x;}
        public void E(double a) {x = a;}
        public void E(double[] a) {x = a[0];}
        public void E(int[] a) {x = a[0];}
        public void E(int i, double a) {x = a;}  //assumes i = 0
        public void Ea1Tv1(double a1, etomica.space.Vector u) {Vector1D u1=(Vector1D)u; x = a1*u1.x;}
        public void Ev1Pa1Tv2(Vector v1, double a1, Vector v2) {
            x = ((Vector1D)v1).x + a1*((Vector1D)v2).x; 
        }
        public void PEa1Tv1(double a1, etomica.space.Vector u) {Vector1D u1=(Vector1D)u; x += a1*u1.x;}
        public void PE(Vector1D u) {x += u.x;}
        public void PE(double a) {x += a;}
        public void ME(Vector1D u) {x -= u.x;}
        public void PE(int i, double a) {x += a;}
        public void TE(double a) {x *= a;}
        public void TE(Vector1D u) {x *= u.x;}
        public void TE(int i, double a) {x *= a;}
        public void DE(double a) {x /= a;}
        public void DE(Vector1D u) {x /= u.x;}
        public void Ev1Pv2(etomica.space.Vector u1, etomica.space.Vector u2) {
            Vector1D v1 = (Vector1D)u1; Vector1D v2 = (Vector1D)u2;
            x = v1.x + v2.x;
        }
        public void Ev1Mv2(etomica.space.Vector u1, etomica.space.Vector u2) {
            Vector1D v1 = (Vector1D)u1; Vector1D v2 = (Vector1D)u2;
            x = v1.x - v2.x;
        }
		public double Mv1Squared(etomica.space.Vector u1) {
			double dx = x-((Vector1D)u1).x;
			return dx*dx;
		}
        public void mod(etomica.space.Vector u) {
            mod((Vector1D)u);
        }
        public void mod(Vector1D u) {
            while(x > u.x) x -= u.x;
            while(x < 0.0) x += u.x;
        }
        public void mod(double a) {
            while(x > a)   x -= a;
            while(x < 0.0) x += a;
        }
//        public void EModShift(Space.Vector r, Space.Vector u) {
//        	EModShift((Vector)r, (Vector)u);
//        }
		//sets this equal to (r mod u) - r
		public void EModShift(Vector r, Vector u) {
            double rx = ((Vector1D)r).x;
            double ux = ((Vector1D)u).x;
            x = rx;
			while(x > ux) x -= ux;
			while(x < 0.) x += ux;
			x -= rx;
		}

        public void EMod2Shift(Vector r, Vector u) {
            double rx = ((Vector1D)r).x;
            double ux = ((Vector1D)u).x;
            x = rx;
            while(x > +ux) x -= (ux+ux);
            while(x < -ux) x += (ux+ux);
            x -= rx;
        }


        public etomica.space.Vector P(Vector u) {Vector1D u1=(Vector1D)u; WORK.x = x+u1.x; return WORK;}
        public etomica.space.Vector M(Vector u) {Vector1D u1=(Vector1D)u; WORK.x = x-u1.x; return WORK;}
        public etomica.space.Vector T(Vector u) {Vector1D u1=(Vector1D)u; WORK.x = x*u1.x; return WORK;}
        public etomica.space.Vector D(Vector u) {Vector1D u1=(Vector1D)u; WORK.x = x/u1.x; return WORK;}
        public void abs() {x = (x>0)?x:-x;}
        public double min() {return x;}
        public double max() {return x;}
        public double squared() {return x*x;}
        public void normalize() {x = 1.0;}
        public double dot(Vector1D u) {return x*u.x;}
        public void transform(etomica.space.Tensor A) {transform((Tensor1D)A);}
        public void transform(Tensor1D A) {x = A.xx * x;}
        public void transform(etomica.space.Boundary b, etomica.space.Vector r0, etomica.space.Tensor A) {transform((Boundary)b, (Vector1D)r0, (Tensor1D)A);}
        public void transform(Boundary b, Vector1D r0, Tensor1D A) {
            WORK.x = x-r0.x; b.nearestImage(WORK); x = r0.x + A.xx*WORK.x;
        }
        public void randomStep(double d) {x += (2.*Simulation.random.nextDouble()-1.0)*d;} //uniformly distributed random step in x and y, within +/- d
        public void setRandom(double d) {x = Simulation.random.nextDouble()*d;}
        public void setRandom(double dx, double dy) {x = Simulation.random.nextDouble()*dx;}
        public void setRandom(Vector1D u) {setRandom(u.x);}
        public void setRandomCube() {x = Simulation.random.nextDouble() - 0.5;}
        public void setRandomInSphere() {setRandomCube();}
        public void setRandomSphere() {randomDirection();}
        public void randomDirection() {x = (Simulation.random.nextDouble() < 0.5) ? -1.0 : +1.0;}
        public void E(etomica.space.Vector u) {E((Vector1D)u);}
        public void PE(etomica.space.Vector u) {PE((Vector1D)u);}
        public void ME(etomica.space.Vector u) {ME((Vector1D)u);}
        public void TE(etomica.space.Vector u) {TE((Vector1D)u);}
        public void DE(etomica.space.Vector u) {DE((Vector1D)u);}
        public double dot(etomica.space.Vector u) {return dot((Vector1D)u);}
        /**
         * Sets this vector equal to its cross product of with a 3D vector.
         * Result is projected into this space, and thus outcome is to make this vector zero.
         */
        public void XE(etomica.space3d.Vector3D u) {x = 0.0;}
        public etomica.space3d.Vector3D cross(etomica.space2d.Vector2D u) {return null;}
        public etomica.space3d.Vector3D cross(etomica.space3d.Vector3D u) {return null;}
        public void randomRotate(double deltheta){//no implementation in 1D
        }
        public double productOfElements() {
            return x;
        }
    }
