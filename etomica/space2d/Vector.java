package etomica.space2d;

import etomica.Simulation;
import etomica.Space2D;
import etomica.Space3D;
import etomica.space.Boundary;
import etomica.space.Tensor;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public final class Vector extends Vector {  //declared final for efficient method calls
        public static final Vector ORIGIN = new Vector(0.0,0.0);  //anything using WORK is not thread-safe
        public static final Vector WORK = new Vector();
        private double x, y;
        public Vector () {x = 0.0; y = 0.0;}
        public Vector (double x, double y) {this.x = x; this.y = y;}
        public Vector (double[] a) {x = a[0]; y = a[1];}//should check length of a for exception
        public Vector (Vector u) {this.E(u);}
        public String toString() {return "("+x+", "+y+")";}
        public double[] toArray() {return new double[] {x, y};}
        public boolean equals(Vector v) {return equals((Vector)v);}
        public boolean equals(Vector v) {return (x == v.x) && (y == v.y);}
        public boolean isZero() {return (x == 0.0) && (y == 0.0);}
        public void sphericalCoordinates(double[] result) {
            result[0] = Math.sqrt(x*x + y*y);
            result[1] = Math.atan2(y,x);  //theta
        }
        public int length() {return Space2D.D;}//bad name for this
        public int D() {return Space2D.D;}
        public double x(int i) {return (i==0) ? x : y;}
        public void setX(int i, double d) {if(i==0) x=d; else y=d;}
        public void E(Vector u) {x = u.x; y = u.y;}
        public void E(double[] u) {x = u[0]; y = u[1];}  //should check length of array for exception
        public void E(double a) {x = a; y = a;}
        public void E(double a, double b) {x = a; y = b;}
        public void Ea1Tv1(double a1, Vector u) {Vector u1=(Vector)u; x = a1*u1.x; y = a1*u1.y;}
        public void PEa1Tv1(double a1, Vector u) {Vector u1=(Vector)u; x += a1*u1.x; y += a1*u1.y;}
        public void PE(Vector u) {x += u.x; y += u.y;}
        public void PE(double a) {x += a; y += a;}
        public void ME(Vector u) {x -= u.x; y -= u.y;}
        public void PE(int i, double a) {if(i==0) x += a; else y += a;}
        public void TE(double a) {x *= a; y *= a;}
        public void TE(Vector u) {x *= u.x; y *= u.y;}
        public void TE(int i, double a) {if(i==0) x *= a; else y *= a;}
        public void DE(double a) {x /= a; y /= a;}
        public void DE(Vector u) {x /= u.x; y /= u.y;}
		public double Mv1Squared(Vector u) {
			Vector u1 = (Vector)u;
			double dx = x-u1.x;
			double dy = y-u1.y;
			return dx*dx + dy*dy;
		}
        public void Ev1Pv2(Vector u1, Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x + v2.x;
            y = v1.y + v2.y;
        }
        public void Ev1Mv2(Vector u1, Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x - v2.x;
            y = v1.y - v2.y;
        }
        public void mod(Vector u) {
            mod((Vector)u);
        }
        public void mod(Vector u) {
            while(x > u.x) x -= u.x;
            while(x < 0.0) x += u.x;
            while(y > u.y) y -= u.y;
            while(y < 0.0) y += u.y;
        }
        public void mod(double a) {
            while(x > a)   x -= a;
            while(x < 0.0) x += a;
            while(y > a)   y -= a;
            while(y < 0.0) y += a;
        }
//		public void EModShift(Space.Vector r, Space.Vector u) {
//			EModShift((Vector)r, (Vector)u);
//		}
		//sets this equal to (r mod u) - r
		public void EModShift(Vector r, Vector u) {
			x = r.x;
			while(x > u.x) x -= u.x;
			while(x < 0.0) x += u.x;
			x -= r.x;
			y = r.y;
			while(y > u.y) y -= u.y;
			while(y < 0.0) y += u.y;
			y -= r.y;
		}


        public Vector P(Vector u) {Vector u1=(Vector)u; WORK.x = x+u1.x; WORK.y = y+u1.y; return WORK;}
        public Vector M(Vector u) {Vector u1=(Vector)u; WORK.x = x-u1.x; WORK.y = y-u1.y; return WORK;}
        public Vector T(Vector u) {Vector u1=(Vector)u; WORK.x = x*u1.x; WORK.y = y*u1.y; return WORK;}
        public Vector D(Vector u) {Vector u1=(Vector)u; WORK.x = x/u1.x; WORK.y = y/u1.y; return WORK;}
        public void abs() {x = (x>0)?x:-x; y = (y>0)?y:-y;}
        public double min() {return (x < y) ? x : y;}
        public double max() {return (x > y) ? x : y;}
        public double squared() {return x*x + y*y;}
        public double dot(Vector u) {return x*u.x + y*u.y;}
        public Vector cross(Vector u) {
        	return new Vector(y*u.x(2), -x*u.x(2), x*u.x(1)-y*u.x(0));
//        	work.setX(0, y*u.x(2));
//            work.setX(1,-x*u.x(2));
//            work.setX(2, x*u.x(1) - y*u.x(0));
        }
        public Vector cross(Vector u) {
			return new Vector(0.0, 0.0, x*u.y - y*u.x);
//            Space3D.Vector.WORK.setX(0, 0.0);
//            Space3D.Vector.WORK.setX(1, 0.0);
//            Space3D.Vector.WORK.setX(2, x*u.y - y*u.x);
//            return Space3D.Vector.WORK;
        }
        /**
         * Replaces this vector with its cross-product with the given 3D vector, with result projected
         * onto the 2D plane.  This vector becomes the result of (this vector) X u.
         */
        public void XE(Vector u) {
            double xNew = y*u.x(2);
            y = -x*u.x(2);
            x = xNew;
        }
            
        public void normalize() {
            double norm = Math.sqrt(1/(x*x + y*y));
            x *= norm;
            y *= norm;
        }
        public void transform(Tensor A) {transform((Tensor)A);}
        public void transform(Tensor A) {
        	double x0 = x;
        	double y0 = y;
            x = A.xx*x0 + A.xy*y0; 
            y = A.yx*x0 + A.yy*y0;
        }
        public void transform(Boundary b, Vector r0, Tensor A) {transform((Boundary)b, (Vector)r0, (Tensor)A);}
        public void transform(Boundary b, Vector r0, Tensor A) {
            WORK.x = x - r0.x; WORK.y = y - r0.y;
            b.nearestImage(WORK);
            x = r0.x + A.xx*WORK.x + A.xy*WORK.y;
            y = r0.y + A.yx*WORK.x + A.yy*WORK.y;
        }
        public void randomStep(double d) {x += (2.*Simulation.random.nextDouble()-1.0)*d; y+= (2.*Simulation.random.nextDouble()-1.0)*d;} //uniformly distributed random step in x and y, within +/- d
        public void setRandom(double d) {x = Simulation.random.nextDouble()*d; y = Simulation.random.nextDouble()*d;}
        public void setRandom(double dx, double dy) {x = Simulation.random.nextDouble()*dx; y = Simulation.random.nextDouble()*dy;}
        public void setRandom(Vector u) {setRandom(u.x,u.y);}
        public void setRandomCube() {
            x = Simulation.random.nextDouble() - 0.5; 
            y = Simulation.random.nextDouble() - 0.5;
        }
        public void setRandomSphere() {
            x = Math.cos(2*Math.PI*Simulation.random.nextDouble()); 
            y = Math.sqrt(1.0 - x*x);
            if(Simulation.random.nextDouble() < 0.5) y = -y;
        }
		// random point in a unit sphere
		public void setRandomInSphere() {//check before using
			double z1 = 0.0;
			double z2 = 0.0;
			double rsq = Double.MAX_VALUE;
			while(rsq > 1.0) {
                
				z1 = 1.0 - 2.0*Simulation.random.nextDouble();
				z2 = 1.0 - 2.0*Simulation.random.nextDouble();
        
				rsq = z1*z1 + z2*z2;
			}
			x = z1;
			y = z2;
		}

        public void randomRotate(double thetaStep){
            double deltheta = (2*Simulation.random.nextDouble() - 1.0)*thetaStep;
            double theta = Math.atan2(y,x);
            theta += deltheta;
            double r = Math.sqrt(x*x + y*y);
            x = r*Math.cos(theta);
            y = r*Math.sin(theta);
        }
        
        public void E(Vector u) {E((Vector)u);}
        public void PE(Vector u) {PE((Vector)u);}
        public void ME(Vector u) {ME((Vector)u);}
        public void TE(Vector u) {TE((Vector)u);}
        public void DE(Vector u) {DE((Vector)u);}
        public double dot(Vector u) {return dot((Vector)u);}
    }