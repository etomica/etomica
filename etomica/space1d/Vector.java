package etomica.space1d;

import etomica.Simulation;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public final class Vector extends etomica.space.Vector {  //declared final for efficient method calls
        public static final Vector ZERO = new Vector(0.0);
        public static final Vector WORK = new Vector();
        double x;
        public Vector () {x = 0.0;}
        public Vector (double a1) {x = a1;}
        public Vector (double[] a) {x = a[0];}//should check length of a for exception
        public Vector (Vector u) {this.E(u);}
        public boolean equals(etomica.space.Vector v) {return equals((Vector)v);}
        public boolean equals(Vector v) {return (x == v.x);}
        public boolean isZero() {return x==0;}
        public String toString() {return "("+x+")";}
        public int length() {return 1;}
        public int D() {return 1;}
        public double x(int i) {return x;}
        public void setX(int i, double d) {x=d;}
        public double[] toArray() {return new double[] {x};}
        public void sphericalCoordinates(double[] result) {result[0] = x;}
        public void E(Vector u) {x = u.x;}
        public void E(double a) {x = a;}
        public void E(double[] a) {x = a[0];}
        public void E(int i, double a) {x = a;}  //assumes i = 0
        public void Ea1Tv1(double a1, etomica.space.Vector u) {Vector u1=(Vector)u; x = a1*u1.x;}
        public void PEa1Tv1(double a1, etomica.space.Vector u) {Vector u1=(Vector)u; x += a1*u1.x;}
        public void PE(Vector u) {x += u.x;}
        public void PE(double a) {x += a;}
        public void ME(Vector u) {x -= u.x;}
        public void PE(int i, double a) {x += a;}
        public void TE(double a) {x *= a;}
        public void TE(Vector u) {x *= u.x;}
        public void TE(int i, double a) {x *= a;}
        public void DE(double a) {x /= a;}
        public void DE(Vector u) {x /= u.x;}
        public void Ev1Pv2(etomica.space.Vector u1, etomica.space.Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x + v2.x;
        }
        public void Ev1Mv2(etomica.space.Vector u1, etomica.space.Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x - v2.x;
        }
		public double Mv1Squared(etomica.space.Vector u1) {
			double dx = x-((Vector)u1).x;
			return dx*dx;
		}
        public void mod(etomica.space.Vector u) {
            mod((Vector)u);
        }
        public void mod(Vector u) {
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
			x = r.x;
			while(x > u.x) x -= u.x;
			while(x < 0.0) x += u.x;
			x -= r.x;
		}


        public etomica.space.Vector P(etomica.space.Vector u) {Vector u1=(Vector)u; WORK.x = x+u1.x; return WORK;}
        public etomica.space.Vector M(etomica.space.Vector u) {Vector u1=(Vector)u; WORK.x = x-u1.x; return WORK;}
        public etomica.space.Vector T(etomica.space.Vector u) {Vector u1=(Vector)u; WORK.x = x*u1.x; return WORK;}
        public etomica.space.Vector D(etomica.space.Vector u) {Vector u1=(Vector)u; WORK.x = x/u1.x; return WORK;}
        public void abs() {x = (x>0)?x:-x;}
        public double min() {return x;}
        public double max() {return x;}
        public double squared() {return x*x;}
        public void normalize() {x = 1.0;}
        public double dot(Vector u) {return x*u.x;}
        public void transform(etomica.space.Tensor A) {transform((Tensor)A);}
        public void transform(Tensor A) {x = A.xx * x;}
        public void transform(etomica.space.Boundary b, etomica.space.Vector r0, etomica.space.Tensor A) {transform((Boundary)b, (Vector)r0, (Tensor)A);}
        public void transform(Boundary b, Vector r0, Tensor A) {
            WORK.x = x-r0.x; b.nearestImage(WORK); x = r0.x + A.xx*WORK.x;
        }
        public void randomStep(double d) {x += (2.*Simulation.random.nextDouble()-1.0)*d;} //uniformly distributed random step in x and y, within +/- d
        public void setRandom(double d) {x = Simulation.random.nextDouble()*d;}
        public void setRandom(double dx, double dy) {x = Simulation.random.nextDouble()*dx;}
        public void setRandom(Vector u) {setRandom(u.x);}
        public void setRandomCube() {x = Simulation.random.nextDouble() - 0.5;}
        public void setRandomInSphere() {setRandomCube();}
        public void setRandomSphere() {randomDirection();}
        public void randomDirection() {x = (Simulation.random.nextDouble() < 0.5) ? -1.0 : +1.0;}
        public void E(etomica.space.Vector u) {E((Vector)u);}
        public void PE(etomica.space.Vector u) {PE((Vector)u);}
        public void ME(etomica.space.Vector u) {ME((Vector)u);}
        public void TE(etomica.space.Vector u) {TE((Vector)u);}
        public void DE(etomica.space.Vector u) {DE((Vector)u);}
        public double dot(etomica.space.Vector u) {return dot((Vector)u);}
        /**
         * Sets this vector equal to its cross product of with a 3D vector.
         * Result is projected into this space, and thus outcome is to make this vector zero.
         */
        public void XE(etomica.space3d.Vector u) {x = 0.0;}
        public etomica.space3d.Vector cross(etomica.space2d.Vector u) {return null;}
        public etomica.space3d.Vector cross(etomica.space3d.Vector u) {return null;}
        public void randomRotate(double deltheta){//no implementation in 1D
        }
    }
