package etomica.space3d;

import etomica.Simulation;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public final class Vector extends etomica.space.Vector {
        public static final Vector ORIGIN = new Vector(0.0, 0.0, 0.0);
        double x, y, z;
        public int length() {return 3;}
        public int D() {return 3;}
        public Vector () {x = 0.0; y = 0.0; z = 0.0;}
        public Vector (double a1, double a2, double a3) {x = a1; y = a2; z = a3;}
        public Vector (double[] a) {x = a[0]; y = a[1]; z = a[2];}//should check length of a for exception
        public Vector (Vector u) {this.E(u);}
        public String toString() {return "("+x+", "+y+", "+z+")";}
        public double x(int i) {return((i==0) ? x : (i==1) ? y : z);}
        public double[] toArray() {return new double[] {x, y, z};}
        public boolean equals(etomica.space.Vector v) {return equals((Vector)v);}
        public boolean equals(Vector v) {return (x == v.x) && (y == v.y) && (z == v.z);}
		public boolean isZero() {return (x == 0.0) && (y == 0.0) && (z == 0);}
        public void sphericalCoordinates(double[] result) {
            result[0] = Math.sqrt(x*x + y*y + z*z);
            result[1] = Math.acos(z/result[0]); //theta
            result[2] = Math.atan2(x, y);  //phi
        }
        public void E(Vector u) {x = u.x; y = u.y; z = u.z;}
        public void E(double a) {x = a; y = a; z = a;}
        public void E(double a, double b, double c) {x = a; y = b; z = c;}
        public void E(int i, double a) {if (i==0) x=a; else if (i==1) y=a; else z=a;}
        public void E(double[] u) {x = u[0]; y = u[1]; z = u[2];}  //should check length of array for exception
        public void Ea1Tv1(double a1, etomica.space.Vector u) {Vector u1=(Vector)u; x = a1*u1.x; y = a1*u1.y; z = a1*u1.z;}
        public void PEa1Tv1(double a1, etomica.space.Vector u) {Vector u1=(Vector)u; x += a1*u1.x; y += a1*u1.y; z += a1*u1.z;}
        public void PE(Vector u) {x += u.x; y += u.y; z += u.z;}
        public void PE(double a) {x += a; y += a; z += a;}
        public void ME(Vector u) {x -= u.x; y -= u.y; z -= u.z;}
        public void PE(int i, double a) {if (i==0) x+=a; else if (i==1) y+=a; else z+=a;}
        public void TE(double a) {x *= a; y *= a; z *= a;}
        public void TE(int i, double a) {if (i==0) x*=a; else if (i==1) y*=a; else z*=a;}
        public void TE(Vector u) {x *= u.x; y *= u.y; z *= u.z;}
        public void DE(double a) {x /= a; y /= a; z /= a;}
        public void DE(Vector u) {x /= u.x; y /= u.y ; z /= u.z;}
        public void Ev1Pv2(etomica.space.Vector u1, etomica.space.Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x + v2.x;
            y = v1.y + v2.y;
            z = v1.z + v2.z;
        }
        public void Ev1Mv2(etomica.space.Vector u1, etomica.space.Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x - v2.x;
            y = v1.y - v2.y;
            z = v1.z - v2.z;
        }
        public void mod(etomica.space.Vector u) {
            mod((Vector)u);
        }
        public void mod(Vector u) {
            while(x > u.x) x -= u.x;
            while(x < 0.0) x += u.x;
            while(y > u.y) y -= u.y;
            while(y < 0.0) y += u.y;
            while(z > u.z) z -= u.z;
            while(z < 0.0) z += u.z;
        }
        public void mod(double a) {
            while(x > a)   x -= a;
            while(x < 0.0) x += a;
            while(y > a)   y -= a;
            while(y < 0.0) y += a;
            while(z > a)   z -= a;
            while(z < 0.0) z += a;
        }
		public void mod2(Vector u) {
			while(x > +u.x) x -= (u.x+u.x);
			while(x < -u.x) x += (u.x+u.x);
			while(y > +u.y) y -= (u.y+u.y);
			while(y < -u.y) y += (u.y+u.y);
			while(z > +u.z) z -= (u.z+u.z);
			while(z < -u.z) z += (u.z+u.z);
		}

//		public void modShift(Vector u, Vector shift) {
//			while(x+shift.x > u.x) shift.x -= u.x;
//			while(x-shift.x < 0.0) shift.x += u.x;
//			while(y+shift.y > u.y) shift.y -= u.y;
//			while(y-shift.y < 0.0) shift.y += u.y;
//			while(z+shift.z > u.z) shift.z -= u.z;
//			while(z-shift.z < 0.0) shift.z += u.z;
//		}
//			don't need Space.Vector form -- this method is used only by methods in Boundary
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
			z = r.z;
			while(z > u.z) z -= u.z;
			while(z < 0.0) z += u.z;
			z -= r.z;
		}
		public void EMod2Shift(Vector r, Vector u) {
			x = r.x;
			while(x > +u.x) x -= (u.x+u.x);
			while(x < -u.x) x += (u.x+u.x);
			x -= r.x;
			y = r.y;
			while(y > +u.y) y -= (u.y+u.y);
			while(y < -u.y) y += (u.y+u.y);
			y -= r.y;
			z = r.z;
			while(z > +u.z) z -= (u.z+u.z);
			while(z < -u.z) z += (u.z+u.z);
			z -= r.z;
		}
        
		public etomica.space.Vector P(etomica.space.Vector u) {etomica.space.Vector work = new Vector(); work.Ev1Pv2(this,u); return work;}
		public etomica.space.Vector M(etomica.space.Vector u) {etomica.space.Vector work = new Vector(); work.Ev1Mv2(this,u); return work;}
		public etomica.space.Vector T(etomica.space.Vector u) {etomica.space.Vector work = new Vector(); work.E(this); work.TE(u); return work;}
		public etomica.space.Vector D(etomica.space.Vector u) {etomica.space.Vector work = new Vector(); work.E(this); work.DE(u); return work;}
		public void abs() {x=(x<0)?-x:x; y=(y<0)?-y:y; z=(z<0)?-z:z;}
        public double min() {return (x < y) ? (x<z)?x:z : (y<z)?y:z;}
        public double max() {return (x > y) ? (x>z)?x:z : (y>z)?y:z;}
        public double squared() {return x*x + y*y + z*z;}
        public double Mv1Pv2Squared(Vector u1, Vector u2) {
        	double dx = x-u1.x+u2.x;
        	double dy = y-u1.y+u2.y;
        	double dz = z-u1.z+u2.z;
        	return dx*dx + dy*dy + dz*dz;
        }
		public double Mv1Squared(etomica.space.Vector u) {
			Vector u1 = (Vector)u;
			double dx = x-u1.x;
			double dy = y-u1.y;
			double dz = z-u1.z;
			return dx*dx + dy*dy + dz*dz;
		}
        public double dot(Vector u) {return x*u.x + y*u.y + z*u.z;}
        public void transform(etomica.space.Tensor A) {transform((Tensor)A);}
        public void transform(Tensor A) {
            double x1 = A.xx*x + A.xy*y + A.xz*z; 
            double y1 = A.yx*x + A.yy*y + A.yz*z;
                    z = A.zx*x + A.zy*y + A.zz*z;
                    x=x1; 
                    y=y1;
        }
        public void transform(etomica.space.Boundary b, etomica.space.Vector r0, etomica.space.Tensor A) {transform((Boundary)b, (Vector)r0, (Tensor)A);}
        public void transform(Boundary b, Vector r0, Tensor A) {
        	this.ME(r0);
        	b.nearestImage(this);
        	this.transform(A);
        	this.PE(r0);
        }
        public void randomStep(double d) {x += (2.*Simulation.random.nextDouble()-1)*d; y +=(2.*Simulation.random.nextDouble()-1)*d; z +=(2.*Simulation.random.nextDouble()-1)*d;}
        public void setRandom(double d) {x = Simulation.random.nextDouble()*d; y = Simulation.random.nextDouble()*d; z = Simulation.random.nextDouble()*d;}
        public void setRandom(double dx, double dy, double dz) {x = Simulation.random.nextDouble()*dx; y = Simulation.random.nextDouble()*dy; z = Simulation.random.nextDouble()*dz;}
        public void setRandom(Vector u) {setRandom(u.x, u.y, u.z);}
        public void setRandomCube() {
            x = Simulation.random.nextDouble() - 0.5;
            y = Simulation.random.nextDouble() - 0.5;
            z = Simulation.random.nextDouble() - 0.5;
        }
        public void setX(int a, double d) { if(a==0) x=d; else if(a==1) y=d; else z=d;}
        //random point on a unit sphere
        public void setRandomSphere() {//check before using
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
        
        // random point in a unit-radius sphere
		public void setRandomInSphere() {//check before using
			double z1 = 0.0;
			double z2 = 0.0;
			double z3 = 0.0;
			double rsq = Double.MAX_VALUE;
			while(rsq > 1.0) {
                
				z1 = 1.0 - 2.0*Simulation.random.nextDouble();
				z2 = 1.0 - 2.0*Simulation.random.nextDouble();
				z3 = 1.0 - 2.0*Simulation.random.nextDouble();
        
				rsq = z1*z1 + z2*z2 + z3*z3;
			}
			x = z1;
			y = z2;
			z = z3;
		}

        /**
        * Creating a random unit vector on unit sphere
        * Uses only two random number generator at a time
        * 
        * @author Jayant Singh
        */
        public void randomVectorOnUnitSphere(){
            double z1=0.0,z2=0.0,zsq=20.0;
            while(zsq > 1.0) {       
                z1  = 2.0*Simulation.random.nextDouble() - 1.0; 
                z2  = 2.0*Simulation.random.nextDouble() - 1.0;  
                zsq = z1*z1 + z2*z2;
            }
                 
            double ranh = 2.0*Math.sqrt(1.0 - zsq) ;   
            x = z1*ranh;     
            y = z2*ranh;     
            z = 1.0 - 2.0*zsq;
        }
        public void randomRotate(double thetaStep){//check before using
            //could be made more efficient by merging with setRandomSphere
            if(thetaStep == 0.0) return;
            if(Math.abs(thetaStep) > Math.PI) {setRandomSphere(); return;}
            double xOld = x; double yOld = y; double zOld = z;
            double r = Math.sqrt(x*x + y*y + z*z);
            double dotMin = r*Math.cos(thetaStep);
            do {
                setRandomSphere();
            } while(xOld*x + yOld*y + zOld*z < dotMin);
            x *= r; y *= r; z*=r;
        }
        
        public void setToOrigin() {x = ORIGIN.x; y = ORIGIN.y; z = ORIGIN.z;}

        public void E(etomica.space.Vector u) {E((Vector) u);}
        public void PE(etomica.space.Vector u) {PE((Vector) u);}
        public void TE(etomica.space.Vector u) {TE((Vector) u);}
        public void ME(etomica.space.Vector u) {ME((Vector) u);}
        public void DE(etomica.space.Vector u) {DE((Vector) u);}
        public double dot(etomica.space.Vector u) {return dot((Vector)u);}
        public Vector cross(etomica.space2d.Vector u) {//does 3d cross-product taking u.z = 0
        	Vector work = new Vector(this);
        	work.x = -z*u.x(1);
            work.y = z*u.x(0);
            work.z = x*u.x(1) - y*u.x(0);
            return work;
        }
        /**
         * Makes a new vector and assigns it the cross product of this with
         * given vector.  This vector is unchanged.
         * @see etomica.space.Vector#cross(Vector)
         */
        public Vector cross(Vector u) {//not thread safe
        	Vector work = new Vector(this);
        	work.XE(u);
        	return work;
        }
        public void XE(Vector u) {//cross product
            double xNew = y*u.z - z*u.y;
            double yNew = z*u.x - x*u.z;
            z = x*u.y - y*u.x;
            y = yNew;
            x = xNew;
        }
        public void normalize() {
            double norm = 1./Math.sqrt(x*x + y*y + z*z);
            x *= norm;
            y *= norm;
            z *= norm;
        }
    }
