package etomica.space2d;

import etomica.Simulation;
import etomica.space.Boundary;
import etomica.space.Vector;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public final class Vector2D extends etomica.space.Vector {  //declared final for efficient method calls
        public static final Vector2D ORIGIN = new Vector2D(0.0,0.0);  //anything using WORK is not thread-safe
        public static final Vector2D WORK = new Vector2D();
        double x, y;
        public Vector2D () {x = 0.0; y = 0.0;}
        public Vector2D (double x, double y) {this.x = x; this.y = y;}
        public Vector2D (double[] a) {x = a[0]; y = a[1];}//should check length of a for exception
        public Vector2D (Vector2D u) {this.E(u);}
        public String toString() {return "("+x+", "+y+")";}
        public void assignTo(double[] array) {array[0] = x; array[1] = y;}
        public double[] toArray() {return new double[] {x, y};}
        public boolean equals(etomica.space.Vector v) {return equals((Vector2D)v);}
        public boolean equals(Vector2D v) {return (x == v.x) && (y == v.y);}
        public boolean isZero() {return (x == 0.0) && (y == 0.0);}
        public void sphericalCoordinates(double[] result) {
            result[0] = Math.sqrt(x*x + y*y);
            result[1] = Math.atan2(y,x);  //theta
        }
        public int length() {return 2;}//bad name for this
        public int D() {return 2;}
        public double x(int i) {return (i==0) ? x : y;}
        public void setX(int i, double d) {if(i==0) x=d; else y=d;}
        public void E(Vector2D u) {x = u.x; y = u.y;}
        public void E(double[] u) {x = u[0]; y = u[1];}  //should check length of array for exception
        public void E(int[] u) {x = u[0]; y = u[1];}  //should check length of array for exception
        public void E(double a) {x = a; y = a;}
        public void E(double a, double b) {x = a; y = b;}
        public void Ea1Tv1(double a1, etomica.space.Vector u) {Vector2D u1=(Vector2D)u; x = a1*u1.x; y = a1*u1.y;}
        public void PEa1Tv1(double a1, etomica.space.Vector u) {Vector2D u1=(Vector2D)u; x += a1*u1.x; y += a1*u1.y;}
        public void PE(Vector2D u) {x += u.x; y += u.y;}
        public void PE(double a) {x += a; y += a;}
        public void ME(Vector2D u) {x -= u.x; y -= u.y;}
        public void PE(int i, double a) {if(i==0) x += a; else y += a;}
        public void TE(double a) {x *= a; y *= a;}
        public void TE(Vector2D u) {x *= u.x; y *= u.y;}
        public void TE(int i, double a) {if(i==0) x *= a; else y *= a;}
        public void DE(double a) {x /= a; y /= a;}
        public void DE(Vector2D u) {x /= u.x; y /= u.y;}
		public double Mv1Squared(etomica.space.Vector u) {
			Vector2D u1 = (Vector2D)u;
			double dx = x-u1.x;
			double dy = y-u1.y;
			return dx*dx + dy*dy;
		}
        public void Ev1Pv2(etomica.space.Vector u1, etomica.space.Vector u2) {
            Vector2D v1 = (Vector2D)u1; Vector2D v2 = (Vector2D)u2;
            x = v1.x + v2.x;
            y = v1.y + v2.y;
        }
        public void Ev1Mv2(etomica.space.Vector u1, etomica.space.Vector u2) {
            Vector2D v1 = (Vector2D)u1; Vector2D v2 = (Vector2D)u2;
            x = v1.x - v2.x;
            y = v1.y - v2.y;
        }
        public void mod(etomica.space.Vector u) {
            mod((Vector2D)u);
        }
        public void mod(Vector2D u) {
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
            Vector2D r2d = (Vector2D)r;
            Vector2D u2d = (Vector2D)u;
			x = r2d.x;
			while(x > u2d.x) x -= u2d.x;
			while(x < 0.0) x += u2d.x;
			x -= r2d.x;
			y = r2d.y;
			while(y > u2d.y) y -= u2d.y;
			while(y < 0.0) y += u2d.y;
			y -= r2d.y;
		}

        public void EMod2Shift(Vector r, Vector u) {
            Vector2D r2d = (Vector2D)r;
            Vector2D u2d = (Vector2D)u;
            x = r2d.x;
            while(x > +u2d.x) x -= (u2d.x+u2d.x);
            while(x < -u2d.x) x += (u2d.x+u2d.x);
            x -= r2d.x;
            y = r2d.y;
            while(y > +u2d.y) y -= (u2d.y+u2d.y);
            while(y < -u2d.y) y += (u2d.y+u2d.y);
            y -= r2d.y;
        }
 
        public etomica.space.Vector P(etomica.space.Vector u) {Vector2D u1=(Vector2D)u; WORK.x = x+u1.x; WORK.y = y+u1.y; return WORK;}
        public etomica.space.Vector M(etomica.space.Vector u) {Vector2D u1=(Vector2D)u; WORK.x = x-u1.x; WORK.y = y-u1.y; return WORK;}
        public etomica.space.Vector T(etomica.space.Vector u) {Vector2D u1=(Vector2D)u; WORK.x = x*u1.x; WORK.y = y*u1.y; return WORK;}
        public etomica.space.Vector D(etomica.space.Vector u) {Vector2D u1=(Vector2D)u; WORK.x = x/u1.x; WORK.y = y/u1.y; return WORK;}
        public void abs() {x = (x>0)?x:-x; y = (y>0)?y:-y;}
        public double min() {return (x < y) ? x : y;}
        public double max() {return (x > y) ? x : y;}
        public double squared() {return x*x + y*y;}
        public double dot(Vector2D u) {return x*u.x + y*u.y;}
        public etomica.space3d.Vector3D cross(etomica.space3d.Vector3D u) {
        	return new etomica.space3d.Vector3D(y*u.x(2), -x*u.x(2), x*u.x(1)-y*u.x(0));
//        	work.setX(0, y*u.x(2));
//            work.setX(1,-x*u.x(2));
//            work.setX(2, x*u.x(1) - y*u.x(0));
        }
        public etomica.space3d.Vector3D cross(Vector2D u) {
			return new etomica.space3d.Vector3D(0.0, 0.0, x*u.y - y*u.x);
//            Space3D.Vector.WORK.setX(0, 0.0);
//            Space3D.Vector.WORK.setX(1, 0.0);
//            Space3D.Vector.WORK.setX(2, x*u.y - y*u.x);
//            return Space3D.Vector.WORK;
        }
        /**
         * Replaces this vector with its cross-product with the given 3D vector, with result projected
         * onto the 2D plane.  This vector becomes the result of (this vector) X u.
         */
        public void XE(etomica.space3d.Vector3D u) {
            double xNew = y*u.x(2);
            y = -x*u.x(2);
            x = xNew;
        }
            
        public void normalize() {
            double norm = Math.sqrt(1/(x*x + y*y));
            x *= norm;
            y *= norm;
        }
        public void transform(etomica.space.Tensor A) {transform((Tensor2D)A);}
        public void transform(Tensor2D A) {
        	double x0 = x;
        	double y0 = y;
            x = A.xx*x0 + A.xy*y0; 
            y = A.yx*x0 + A.yy*y0;
        }
        public void transform(etomica.space.Boundary b, etomica.space.Vector r0, etomica.space.Tensor A) {transform((Boundary)b, (Vector2D)r0, (Tensor2D)A);}
        public void transform(Boundary b, Vector2D r0, Tensor2D A) {
            WORK.x = x - r0.x; WORK.y = y - r0.y;
            b.nearestImage(WORK);
            x = r0.x + A.xx*WORK.x + A.xy*WORK.y;
            y = r0.y + A.yx*WORK.x + A.yy*WORK.y;
        }
        public void randomStep(double d) {x += (2.*Simulation.random.nextDouble()-1.0)*d; y+= (2.*Simulation.random.nextDouble()-1.0)*d;} //uniformly distributed random step in x and y, within +/- d
        public void setRandom(double d) {x = Simulation.random.nextDouble()*d; y = Simulation.random.nextDouble()*d;}
        public void setRandom(double dx, double dy) {x = Simulation.random.nextDouble()*dx; y = Simulation.random.nextDouble()*dy;}
        public void setRandom(Vector2D u) {setRandom(u.x,u.y);}
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
        
        public void E(etomica.space.Vector u) {E((Vector2D)u);}
        public void PE(etomica.space.Vector u) {PE((Vector2D)u);}
        public void ME(etomica.space.Vector u) {ME((Vector2D)u);}
        public void TE(etomica.space.Vector u) {TE((Vector2D)u);}
        public void DE(etomica.space.Vector u) {DE((Vector2D)u);}
        public double dot(etomica.space.Vector u) {return dot((Vector2D)u);}
        public double productOfElements() {
            return x*y;
        }
    }
