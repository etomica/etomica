package etomica;

//centralImage not updated to molecule form as in Space3D

 /* History of changes
  * 09/01/02 (DAK) added accelerateTo method to Coordinate
  *                changed CoordinateGroup.randomizeMomentum to not enforce zero COM momentum
  * 09/05/02 (DAK) fixed error in accelerateTo (still probably does not do what one expects
  *                if accelerating to nonzero momentum).
  * 01/04/03 (SKK/DAK) added HSlit Boundary class for horizontal-slit PBC.
  * 01/12/03 (JKS/DAK) corrected error in Vector.transform, where updated xyz
  * 07/10/03 (DAK) added resetV method to CoordinatePair
  * 08/27/03 (DAK) added isZero method to Vector
  * 08/29/03 (DAK) implemented centralImage(Space.Coordinate) in Boundary
  * 12/09/03 (DAK) added setRandomInSphere method in Vector
  * 01/22/04 (DAK) added positionCOM and translateCOMTo to CoordinateGroup;
  * redefined position() in CoordinateGroup to be first-atom position (as it has
  * been for Space3D for some time now).
  */
public class Space2D extends Space implements EtomicaElement {
    
    public static final int D = 2;
    public final int D() {return D;}
    public final int powerD(int n) {return n*n;}
    public final double powerD(double a) {return a*a;}
    public static final Vector ORIGIN = new Vector();
    public final Space.Vector origin() {return ORIGIN;}
    public static final Space2D INSTANCE = new Space2D();
    
    public Space2D() {super(2);}
    
    public double sphereVolume(double r) {return Math.PI*r*r;}  //volume of a sphere of radius r
    public double sphereArea(double r) {return 2.0*Math.PI*r;}  //surface area of sphere of radius r (used for differential shell volume)
    public Space.Vector makeVector() {return new Vector();}
    public Space.Orientation makeOrientation() {return new Orientation();}
    public Space.Tensor makeTensor() {return new Tensor();}
    public Space.Tensor makeRotationTensor() {return new RotationTensor();}
    public Space.Coordinate makeCoordinate(Atom a) {
        if(a.node instanceof AtomTreeNodeGroup) return new CoordinateGroup(a);
        else if(a.type instanceof AtomType.Rotator) return new OrientedCoordinate(a);
        else {return new Coordinate(a);}
    }
    public Space.CoordinatePair makeCoordinatePair() {return new CoordinatePair();}
    
    public Space.Boundary.Type[] boundaryTypes() {return Boundary.TYPES;}
    public Space.Boundary makeBoundary() {return makeBoundary(Boundary.PERIODIC_SQUARE);}  //default
    public Space.Boundary makeBoundary(Space.Boundary.Type t) {
        if(t == Boundary.NONE) {return new BoundaryNone();}
        else if(t == Boundary.PERIODIC_SQUARE) {return new BoundaryPeriodicSquare();}
        else if(t == Boundary.HSLIT) {return new BoundaryHSlit();}
 //       else if(t == Boundary.SLIDING_BRICK) return new BoundarySlidingBrick();
        else return null;
    }
    
    public int[] makeArrayD(int i) {return new int[] {i, i};}
    public double[] makeArrayD(double d) {return new double[] {d, d};}
 
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Two-dimensional space");
        return info;
    }

    public static final double r2(Vector u1, Vector u2, Boundary b) {
        Vector.WORK.x = u1.x - u2.x;
        Vector.WORK.y = u1.y - u2.y;
        b.nearestImage(Vector.WORK);
        return Vector.WORK.x*Vector.WORK.x + Vector.WORK.y*Vector.WORK.y;
    }
        
    public final static class Vector extends Space.Vector {  //declared final for efficient method calls
        public static final Vector ORIGIN = new Vector(0.0,0.0);  //anything using WORK is not thread-safe
        public static final Vector WORK = new Vector();
        private double x, y;
        public Vector () {x = 0.0; y = 0.0;}
        public Vector (double x, double y) {this.x = x; this.y = y;}
        public Vector (double[] a) {x = a[0]; y = a[1];}//should check length of a for exception
        public Vector (Vector u) {this.E(u);}
        public String toString() {return "("+x+", "+y+")";}
        public double[] toArray() {return new double[] {x, y};}
        public boolean equals(Space.Vector v) {return equals((Vector)v);}
        public boolean equals(Vector v) {return (x == v.x) && (y == v.y);}
        public boolean isZero() {return (x == 0.0) && (y == 0.0);}
        public void sphericalCoordinates(double[] result) {
            result[0] = Math.sqrt(x*x + y*y);
            result[1] = Math.atan2(y,x);  //theta
        }
        public int length() {return D;}//bad name for this
        public int D() {return D;}
        public double x(int i) {return (i==0) ? x : y;}
        public void setX(int i, double d) {if(i==0) x=d; else y=d;}
        public void E(Vector u) {x = u.x; y = u.y;}
        public void E(double[] u) {x = u[0]; y = u[1];}  //should check length of array for exception
        public void E(double a) {x = a; y = a;}
        public void E(double a, double b) {x = a; y = b;}
        public void Ea1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x = a1*u1.x; y = a1*u1.y;}
        public void PEa1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x += a1*u1.x; y += a1*u1.y;}
        public void PE(Vector u) {x += u.x; y += u.y;}
        public void PE(double a) {x += a; y += a;}
        public void ME(Vector u) {x -= u.x; y -= u.y;}
        public void PE(int i, double a) {if(i==0) x += a; else y += a;}
        public void TE(double a) {x *= a; y *= a;}
        public void TE(Vector u) {x *= u.x; y *= u.y;}
        public void TE(int i, double a) {if(i==0) x *= a; else y *= a;}
        public void DE(double a) {x /= a; y /= a;}
        public void DE(Vector u) {x /= u.x; y /= u.y;}
		public double Mv1Squared(Space.Vector u) {
			Vector u1 = (Vector)u;
			double dx = x-u1.x;
			double dy = y-u1.y;
			return dx*dx + dy*dy;
		}
        public void Ev1Pv2(Space.Vector u1, Space.Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x + v2.x;
            y = v1.y + v2.y;
        }
        public void Ev1Mv2(Space.Vector u1, Space.Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x - v2.x;
            y = v1.y - v2.y;
        }
        public void mod(Space.Vector u) {
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


        public Space.Vector P(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x+u1.x; WORK.y = y+u1.y; return WORK;}
        public Space.Vector M(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x-u1.x; WORK.y = y-u1.y; return WORK;}
        public Space.Vector T(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x*u1.x; WORK.y = y*u1.y; return WORK;}
        public Space.Vector D(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x/u1.x; WORK.y = y/u1.y; return WORK;}
        public void abs() {x = (x>0)?x:-x; y = (y>0)?y:-y;}
        public double min() {return (x < y) ? x : y;}
        public double max() {return (x > y) ? x : y;}
        public double squared() {return x*x + y*y;}
        public double dot(Vector u) {return x*u.x + y*u.y;}
        public Space3D.Vector cross(Space3D.Vector u) {
        	return new Space3D.Vector(y*u.x(2), -x*u.x(2), x*u.x(1)-y*u.x(0));
//        	work.setX(0, y*u.x(2));
//            work.setX(1,-x*u.x(2));
//            work.setX(2, x*u.x(1) - y*u.x(0));
        }
        public Space3D.Vector cross(Space2D.Vector u) {
			return new Space3D.Vector(0.0, 0.0, x*u.y - y*u.x);
//            Space3D.Vector.WORK.setX(0, 0.0);
//            Space3D.Vector.WORK.setX(1, 0.0);
//            Space3D.Vector.WORK.setX(2, x*u.y - y*u.x);
//            return Space3D.Vector.WORK;
        }
        /**
         * Replaces this vector with its cross-product with the given 3D vector, with result projected
         * onto the 2D plane.  This vector becomes the result of (this vector) X u.
         */
        public void XE(Space3D.Vector u) {
            double xNew = y*u.x(2);
            y = -x*u.x(2);
            x = xNew;
        }
            
        public void normalize() {
            double norm = Math.sqrt(1/(x*x + y*y));
            x *= norm;
            y *= norm;
        }
        public void transform(Space.Tensor A) {transform((Tensor)A);}
        public void transform(Tensor A) {
        	double x0 = x;
        	double y0 = y;
            x = A.xx*x0 + A.xy*y0; 
            y = A.yx*x0 + A.yy*y0;
        }
        public void transform(Space.Boundary b, Space.Vector r0, Space.Tensor A) {transform((Boundary)b, (Vector)r0, (Tensor)A);}
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
        
        public void E(Space.Vector u) {E((Vector)u);}
        public void PE(Space.Vector u) {PE((Vector)u);}
        public void ME(Space.Vector u) {ME((Vector)u);}
        public void TE(Space.Vector u) {TE((Vector)u);}
        public void DE(Space.Vector u) {DE((Vector)u);}
        public double dot(Space.Vector u) {return dot((Vector)u);}
    }
    
    public static class Tensor implements Space.Tensor {
        public int length() {return D;}
        double xx, xy, yx, yy;
        public static final Tensor ZERO = new Tensor();
        public static final Tensor WORK = new Tensor();  //anything using WORK is not thread-safe
        public Tensor () {xx = xy = yx = yy = 0.0;}
        public Tensor (double[] d) {
            this.E(d);
        }

        public double component(int i, int j) {
            return (i==0) ? ( (j==0) ? xx : xy ) : ( (j==0) ? yx : yy );}
        public void setComponent(int i, int j, double d) {
            if(i==0) {if(j==0) xx=d; else xy=d;}
            else     {if(j==0) yx=d; else yy=d;}
        }
        public void E(Tensor t) {xx=t.xx; xy=t.xy; yx=t.yx; yy=t.yy;}
        public void E(Vector u1, Vector u2) {xx=u1.x*u2.x; xy=u1.x*u2.y; yx=u1.y*u2.x; yy=u1.y*u2.y;}
        public void E(double a) {xx = xy = yx = yy = a;}
        public void PE(Tensor t) {xx+=t.xx; xy+=t.xy; yx+=t.yx; yy+=t.yy;}
        public void PE(int i, int j, double a) {
            if(i==0) {if(j==0) xx+=a; else xy+=a;}
            else     {if(j==0) yx+=a; else yy+=a;}
        }
        public void PE(Vector u1, Vector u2) {xx+=u1.x*u2.x; xy+=u1.x*u2.y; yx+=u1.y*u2.x; yy+=u1.y*u2.y;}
        public double trace() {return xx + yy;}
        
        public void E(Space.Tensor t) {E((Tensor)t);}
        public void E(Space.Vector u1, Space.Vector u2) {E((Vector)u1, (Vector)u2);}
        public void PE(Space.Tensor t) {PE((Tensor)t);}
        public void PE(Space.Vector u1, Space.Vector u2) {PE((Vector)u1, (Vector)u2);}
        public void TE(double a) {xx*=a; xy*=a; yx*=a; yy*=a;}
        public void E(double[] d) {
            if(d.length != 4) throw new IllegalArgumentException("Array size incorrector for tensor");
            xx = d[0]; xy = d[1]; 
            yx = d[2]; yy = d[3];
        }
        public void assignTo(double[] d) {
            if(d.length != 1) throw new IllegalArgumentException("Array size incorrector for tensor");
            d[0] = xx; d[1] = xy; 
            d[2] = yx; d[3] = yy;
        }

    }
    
    public static class RotationTensor extends Tensor implements Space.RotationTensor {
        public RotationTensor() {super(); reset();}
        public void reset() {
            xx = 1.0; xy = 0.0;
            yx = 0.0; yy = 1.0;
        }
        public void setAxial(int i, double theta) {
            double st = Math.sin(theta);
            double ct = Math.cos(theta);
            switch(i) {
                case 2: xx = ct; xy=-st;
                        yx = st; yy=ct;
                        return;
                default: throw new IllegalArgumentException();
            }
        }
        public void setAngles(double[] angles) {}
        public void invert() {xy *= -1; yx *= -1;}
    }

    public static final class CoordinatePair extends Space.CoordinatePair {  
        Coordinate c1;
        Coordinate c2;
        private final Vector dr = new Vector();  //note that dr is not cloned if this is cloned -- should change this if using dr in clones; also this makes cloned coordinatePairs not thread-safe
        private double drx, dry, dvx, dvy;
		private Space.Boundary boundary = Space.Boundary.NULL;
        public CoordinatePair() {super();}
        public double r2() {return r2;}
		public void setBoundary(Space.Boundary b) {this.boundary = b;}
		public Space.Boundary getBoundary() {return boundary;}		
        public void reset(Space.Coordinate coord1, Space.Coordinate coord2) {  //don't usually use this; instead set c1 and c2 directly, without a cast
            c1 = (Coordinate)coord1;
            c2 = (Coordinate)coord2;
            reset();
        }
        public void reset() {
            dr.x = c2.r.x - c1.r.x;
            dr.y = c2.r.y - c1.r.y;
            boundary.nearestImage(dr);
//            c1.atom.node.parentPhase().boundary().nearestImage(dr);
            drx = dr.x; 
            dry = dr.y;
            r2 = drx*drx + dry*dry;
        }
        public void resetV() {
            double rm1 = c1.rm();
            double rm2 = c2.rm();
            dvx = (rm2*c2.p.x - rm1*c1.p.x);  
            dvy = (rm2*c2.p.y - rm1*c1.p.y);  
        }
        /**
         * Recomputes pair separation, with atom 2 shifted by the given vector
         * Does not apply any PBC, regardless of boundary chosen for space
         */
        public void reset(Space2D.Vector M) {
            dr.x = c2.r.x - c1.r.x + M.x;
            dr.y = c2.r.y - c1.r.y + M.y;
            drx = dr.x;
            dry = dr.y;
            r2 = drx*drx + dry*dry;
        }
        public Space.Vector dr() {return dr;}
        public double dr(int i) {return (i==0) ? drx : dry;}
        public double dv(int i) {return (i==0) ? dvx : dvy;}
        public double v2() {
            return dvx*dvx + dvy*dvy;
        }
        public double vDot(Space.Vector u) {return vDot((Space2D.Vector)u);}
        public double vDot(Space2D.Vector u) {return dvx*u.x + dvy*u.y;}
        public double vDotr() {
            return drx*dvx + dry*dvy;
        }
        public void push(double impulse) {  //changes momentum in the direction joining the atoms
            c1.p.x += impulse*drx;
            c1.p.y += impulse*dry;
            c2.p.x -= impulse*drx;
            c2.p.y -= impulse*dry;
        }
        public void nudge(double rDelta) {
            double ratio = c2.mass()*c1.rm()*rDelta;
            c1.r.x -= ratio*dr.x;
            c1.r.y -= ratio*dr.y;
            c2.r.x += ratio*dr.x;
            c2.r.y += ratio*dr.y;
        }
        public void setSeparation(double r2New) {
            double ratio = c2.mass()*c1.rm();  // (mass2/mass1)
            double delta = (Math.sqrt(r2New/this.r2()) - 1.0)/(1+ratio);
            c1.r.x -= ratio*delta*drx;
            c1.r.y -= ratio*delta*dry;
            c2.r.x += delta*drx;
            c2.r.y += delta*dry;
            //need call reset?
        }
    }

    public static class Coordinate extends Space.Coordinate {
        public Coordinate nextCoordinate, previousCoordinate;
        public final Vector r = new Vector();  //Cartesian coordinates
        public final Vector p = new Vector();  //Momentum vector
        public final Vector rLast = new Vector();  //vector for saving position
        protected final Vector work = new Vector();
        public Coordinate(Atom a) {super(a);}
        
        public void setNextAtom(Atom a) {
            if(a == null) nextCoordinate = null;
            else {
                nextCoordinate = (Coordinate)a.coord;
                ((Coordinate)a.coord).previousCoordinate = this;
            }
        }
        public Atom nextAtom() {return nextCoordinate!=null ? nextCoordinate.atom : null;}
        public Atom previousAtom() {return previousCoordinate!=null ? previousCoordinate.atom : null;}
        public void clearPreviousAtom() {previousCoordinate = null;}
                
        public Space.Vector position() {return r;}
        public Space.Vector momentum() {return p;}
        public double position(int i) {return r.x(i);}
        public double momentum(int i) {return p.x(i);}
        public double kineticEnergy() {return 0.5*p.squared()*rm();}
        public void freeFlight(double t) {
            double tM = t*rm(); // t/mass
            r.x += p.x*tM;
            r.y += p.y*tM;
        }
        
        public void inflate(double s) {r.x *= s; r.y *= s;}
        public void inflate(Space.Vector s) {Vector u = (Vector)s; r.x *= u.x; r.y *= u.y;}
        
        public void transform(Space.Vector r0, Space.Tensor A) {
            r.transform((Boundary)atom.node.parentPhase().boundary(), (Vector)r0, (Tensor)A);
        }
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateBy(Space.Vector u) {
            r.PE((Space2D.Vector)u);
        }
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateBy(double d, Space.Vector u) {
            r.PEa1Tv1(d,(Vector)u);
        }
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateTo(Space.Vector u) {
            r.E((Space2D.Vector)u);
        }      
        public void displaceBy(Space.Vector u) {rLast.E(r); translateBy((Vector)u);}
        public void displaceBy(double d, Space.Vector u) {rLast.E(r); translateBy(d,(Vector)u);}
        public void displaceTo(Space.Vector u) {rLast.E(r); translateTo((Vector)u);}  
        public void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        public void displaceToRandom(etomica.Phase p) {rLast.E(r); translateToRandom(p);}
        public void replace() {
            r.E(rLast);
        }
    //    public final void inflate(double s) {r.TE(s);}

        public void accelerateBy(Space.Vector u) {p.PE(u);}
        public void accelerateBy(double d, Space.Vector u) {p.PEa1Tv1(d,u);}
        public void accelerateTo(Space.Vector u) {p.E(u);}

        public void randomizeMomentum(double temperature) {
            if(isStationary()) {p.E(0.0); return;}
            p.setX(0,MaxwellBoltzmann.randomMomentumComponent(temperature,mass()));
            p.setX(1,MaxwellBoltzmann.randomMomentumComponent(temperature,mass()));
        }
    }
    
    public static class CoordinateGroup extends Coordinate {
        
        private final AtomIteratorListSimple childIterator = new AtomIteratorListSimple();
		private Atom firstAtom;
        
        public CoordinateGroup(Atom a) {
            super(a);
            childIterator.setList(((AtomTreeNodeGroup)a.node).childList);
        }

        /**
         * Applies transformation to COM of group, keeping all internal atoms at same relative
         * positions.
         */
        public void transform(Space.Vector r0, Space.Tensor A) {
            work.E(position()); //work = r
            work.transform(atom.node.parentPhase().boundary(), r0, A);
            work.ME(r);//now work vector contains translation vector for COM
            translateBy(work);
        }
        public Space.Vector position() {
			if(firstAtom == null) firstAtom = ((AtomTreeNodeGroup)atom.node).childList.getFirst(); //DAK 
			if(firstAtom == null) {r.E(0.0); return r;}
			return firstAtom.coord.position();
		}
		public Space.Vector positionCOM() {
            r.E(0.0); double massSum = 0.0;
            childIterator.reset();
            while(childIterator.hasNext()) {
                Atom a = childIterator.nextAtom();
                r.PEa1Tv1(a.coord.mass(), a.coord.positionCOM()); 
                massSum += a.coord.mass();
            }
            r.DE(massSum);
            return r;
        }
        public Space.Vector momentum() {
            p.E(0.0);
            childIterator.reset();
            while(childIterator.hasNext()) {
                p.PE(childIterator.nextAtom().coord.momentum());
            }
            return p;
        }
        public double position(int i) {
            double sum = 0.0; double massSum = 0.0;
            childIterator.reset();
            while(childIterator.hasNext()) {
                Atom a = childIterator.nextAtom();
                sum += a.coord.mass()*a.coord.position(i); 
                massSum += a.coord.mass();
            }
            sum /= massSum;
            return sum;
        }
        public double momentum(int i) {
            double sum = 0.0;
            childIterator.reset();
            while(childIterator.hasNext()) {
                Atom a = childIterator.nextAtom();
                sum += a.coord.mass()*a.coord.momentum(i);
            }
            return sum;
        }
        public double kineticEnergy() {
            double sum = 0.0;
            childIterator.reset();
            while(childIterator.hasNext()) {
                sum += childIterator.nextAtom().coord.kineticEnergy();
            }
            return sum;
        }
        public void freeFlight(double t) {
            double sum = 0.0;
            childIterator.reset();
            while(childIterator.hasNext()) {
                childIterator.nextAtom().coord.freeFlight(t);
            }
        }
        public void inflate(double scale) {
            work.E(position());
            work.TE(scale-1.0);
            translateBy(work);
        }
        public void inflate(Space.Vector scale) {
            scale.PE(-1.0);
            work.E(position());
            work.TE(scale);
            translateBy(work);
            scale.PE(1.0);
        }
        public void translateBy(Space.Vector u) {
            childIterator.reset();
            while(childIterator.hasNext()) {
                childIterator.nextAtom().coord.translateBy(u);
            }
        }
        public void translateBy(double d, Space.Vector u) {
            childIterator.reset();
            while(childIterator.hasNext()) {
                childIterator.nextAtom().coord.translateBy(d, u);
            }
        }
		/**
		 * Translates center of mass of group to the given position.
		 * @param u new position for the COM
		 */
		public void translateCOMTo(Space.Vector u) {
			work.Ea1Tv1(-1,positionCOM()); //position() uses work, so need this first
			work.PE(u);
			translateBy(work);
		}
		/**
		 * Translates group so that first atom is at the given position.
		 */
        public void translateTo(Space.Vector u) {
            work.Ea1Tv1(-1,position());
            work.PE(u);
            translateBy(work);
        }
        //plan to do away with displace methods
        public void displaceBy(Space.Vector u) {
            childIterator.reset();
            while(childIterator.hasNext()) {
                Atom a = childIterator.nextAtom();
                a.coord.displaceBy(u);
            }
        }
        public void displaceBy(double d, Space.Vector u) {
            childIterator.reset();
            while(childIterator.hasNext()) {
                Atom a = childIterator.nextAtom();
                a.coord.displaceBy(d, u);
            }
        }
        public void displaceTo(Space.Vector u) {
            work.Ea1Tv1(-1,position()); //position() uses work, so need this first
            work.PE(u);
            displaceBy(work);
        }
        public void displaceToRandom(etomica.Phase p) {
            displaceTo(p.boundary().randomPosition());
        }
        public void replace() {
            childIterator.reset();
            while(childIterator.hasNext()) {
                Atom a = childIterator.nextAtom();
                a.coord.replace();
            }
        }
        public void accelerateBy(Space.Vector u) {
            childIterator.reset();
            while(childIterator.hasNext()) {
                Atom a = childIterator.nextAtom();
                a.coord.accelerateBy(u);
            }
        }
        public void accelerateBy(double d, Space.Vector u) {
            childIterator.reset();
            while(childIterator.hasNext()) {
                Atom a = childIterator.nextAtom();
                a.coord.accelerateBy(d, u);
            }
        }
        public void accelerateTo(Space.Vector u) {
            work.Ea1Tv1(-1.0/childIterator.size(),momentum());//probably need this first
            work.PE(u);
            accelerateBy(work);
        }
        
        public final void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
            
        public void randomizeMomentum(double temperature) {
            switch(((AtomTreeNodeGroup)atom.node).childAtomCount()) {
                case 0: return;
                case 1: ((AtomTreeNodeGroup)atom.node).firstChildAtom().coord.randomizeMomentum(temperature);//do not zero COM momentum if only one child atom
                        return;
                default://multi-atom group
     //               work.E(0.0); double sum=0.0;
                    childIterator.reset();
                    while(childIterator.hasNext()) {
                        Atom a = childIterator.nextAtom();
                        a.coord.randomizeMomentum(temperature);
     //                   work.PE(a.coord.momentum());
     //                   sum++;
                    }
     //               work.DE(-sum);
     //               childIterator.reset();
     //               while(childIterator.hasNext()) {
     //                   Atom a = childIterator.next();
   //                     a.coord.accelerateBy(work);
     //               }
            }//end switch
        }//end randomizeMomentum
            
    }//end of CoordinateGroup

    public static class OrientedCoordinate extends Coordinate implements Space.Coordinate.Angular {
        private double L = 0.0; //magnitude of angular momentum
        private final Space3D.Vector vector = new Space3D.Vector();//used to return vector quantities (be sure to keep x and y components zero)
        private final double[] I;
        private final Orientation orientation = new Orientation();
        public OrientedCoordinate(Atom a) {
            super(a);
            I = ((AtomType.SphericalTop)((Atom)a).type).momentOfInertia();
        }
        public Space3D.Vector angularMomentum() {vector.setX(2,L); return vector;}
        public Space3D.Vector angularVelocity() {vector.setX(2,L/I[0]); return vector;}
        public void angularAccelerateBy(Space3D.Vector t) {L += t.x(2);}
        public Space.Orientation orientation() {return orientation;}
        public double kineticEnergy() {return super.kineticEnergy() + 0.5*L*L/I[0];}
        public void freeFlight(double t) {
            super.freeFlight(t);
            orientation.rotateBy(t*L/I[0]);//all elements of I equal for spherical top
        }
    }
    public static class Orientation extends Space.Orientation {
        //The rotation matrix A operates on the components of a vector in the space-fixed frame to yield the
        //components in the body-fixed frame
        private final double[][] A = new double[D][D];
        private final Vector[] bodyFrame = new Vector[] {new Vector(1.0,0.0), new Vector(0.0,1.0)};
        private final double[] angle = new double[1];
        private boolean needToUpdateA = true;
        public void E(Space.Orientation o) {E((Orientation)o);}
        public void E(Orientation o) {angle[0] = o.angle[0]; needToUpdateA = true;}
        public Space.Vector[] bodyFrame() {
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
        public void convertToBodyFrame(Space.Vector v) {convertToBodyFrame((Vector)v);}
        public void convertToSpaceFrame(Space.Vector v) {convertToSpaceFrame((Vector)v);}
    }
    
    public static abstract class Boundary extends Space.Boundary {
        public static class Type extends Space.Boundary.Type {
            public Type(String label) {super(label);}
            public Constants.TypedConstant[] choices() {return TYPES;}
        }
        public static final Type NONE = new Type("None");
        public static final Type PERIODIC_SQUARE = new Type("Periodic Square");
        public static final Type HSLIT = new Type("Periodic HSlit");
//        public static final Type SLIDING_BRICK = new Type("Sliding Brick");
        public static final Type[] TYPES = {NONE,PERIODIC_SQUARE, HSLIT/*,SLIDING_BRICK*/};
        public Boundary() {super();}
        public Boundary(Phase p) {super(p);}
        public abstract void nearestImage(Vector dr);
        public abstract boolean centralImage(Vector r);
        public abstract boolean centralImage(Coordinate c);
        public boolean centralImage(Space.Coordinate c) {return centralImage((Coordinate)c);}
    }

    /**
     * Class for implementing no periodic boundary conditions
     */
    protected static class BoundaryNone extends Boundary {
        private final Vector temp = new Vector();
        private final Vector dimensions = new Vector(Default.BOX_SIZE, Default.BOX_SIZE);
        private final Vector dimensionsCopy = new Vector();
        public Space.Vector dimensions() {dimensionsCopy.E(dimensions); return dimensionsCopy;}
        public BoundaryNone() {super();}
        public BoundaryNone(Phase p) {super(p);}
        public Space.Boundary.Type type() {return Boundary.NONE;}
        public final void nearestImage(Space.Vector dr) {}
        public final boolean centralImage(Space.Vector r) {return false;}
        public final void nearestImage(Vector dr) {}
        public final boolean centralImage(Vector r) {return false;}
        public final boolean centralImage(Coordinate c) {return false;}
        public double volume() {return dimensions.x*dimensions.y;}
        public void inflate(double s) {dimensions.TE(s);}
        public void inflate(Space.Vector s) {dimensions.TE(s);}
        public void setDimensions(Space.Vector v) {dimensions.E(v);}
        public double[][] imageOrigins(int nShells) {return new double[0][D];}
        public float[][] getOverflowShifts(Space.Vector rr, double distance) {return shift0;}
        public Space.Vector randomPosition() {  //arbitrary choice for this method in this boundary
            temp.x = dimensions.x*Simulation.random.nextDouble(); 
            temp.y = dimensions.y*Simulation.random.nextDouble(); 
            return temp;
        }
    }//end of BoundaryNone

    /**
     * Class for implementing rectangular periodic boundary conditions
     */
    protected static class BoundaryPeriodicSquare extends Boundary implements Space.Boundary.Periodic {
        public BoundaryPeriodicSquare() {this(Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p) {this(p,Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p, double lx, double ly) {super(p); dimensions.x = lx; dimensions.y = ly; updateDimensions();}
        public BoundaryPeriodicSquare(double lx, double ly) {dimensions.x = lx; dimensions.y = ly; updateDimensions();}
        public Space.Boundary.Type type() {return Boundary.PERIODIC_SQUARE;}
        private final Vector temp = new Vector();
        private final Vector dimensions = new Vector();
        private final Vector dimensionsCopy = new Vector();
        private final Vector dimensionsHalf = new Vector();
        public final Space.Vector dimensions() {return dimensionsCopy;}
        //temporary methods included to test combobox selector
        public Constants.Alignment align = Constants.VERTICAL;
        public Constants.Alignment getAlign() {return align;}
        public void setAlign(Constants.Alignment newAlign) {align = newAlign;}
        //end of temporary methods
        private final void updateDimensions() {
            dimensionsHalf.Ea1Tv1(0.5,dimensions);
            dimensionsCopy.E(dimensions);
        }
        public Space.Vector randomPosition() {
            temp.x = dimensions.x*Simulation.random.nextDouble(); 
            temp.y = dimensions.y*Simulation.random.nextDouble(); 
            return temp;
        }
        public void nearestImage(Space.Vector dr) {nearestImage((Vector)dr);}
        public void nearestImage(Vector dr) {
           // dr.x -= dimensions.x * ((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x+0.5) : Math.ceil(dr.x/dimensions.x-0.5));
           // dr.y -= dimensions.y * ((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y+0.5) : Math.ceil(dr.y/dimensions.y-0.5));
           // final double dimxHalf = 0.5*dimensions.x;
           // final double dimyHalf = 0.5*dimensions.y;
            while(dr.x > +dimensionsHalf.x) dr.x -= dimensions.x;
            while(dr.x < -dimensionsHalf.x) dr.x += dimensions.x;
            while(dr.y > +dimensionsHalf.y) dr.y -= dimensions.y;
            while(dr.y < -dimensionsHalf.y) dr.y += dimensions.y;
        }
        public boolean centralImage(Coordinate c) {return centralImage(c.r);}
        public boolean centralImage(Space.Vector r) {return centralImage((Vector)r);}
        public boolean centralImage(Vector r) {
        /*    while(r.x > dimensions.x) r.x -= dimensions.x;
            while(r.x < 0.0)          r.x += dimensions.x;
            while(r.y > dimensions.y) r.y -= dimensions.y;
            while(r.y < 0.0)          r.y += dimensions.y;
        */  boolean changed = false;
            while(r.x > dimensions.x) {r.x -= dimensions.x; changed = true;}
            while(r.x < 0.0)          {r.x += dimensions.x; changed = true;}
            while(r.y > dimensions.y) {r.y -= dimensions.y; changed = true;}
            while(r.y < 0.0)          {r.y += dimensions.y; changed = true;}
            return changed;  
            //r.x -= dimensions.x * ((r.x >= 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
            //r.y -= dimensions.y * ((r.y >= 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y-1.0));
        }
        public void inflate(double scale) {
            dimensions.TE(scale); 
            updateDimensions();
            phase().boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
        }
        public void inflate(Space.Vector scale) {
            dimensions.TE(scale); 
            updateDimensions();
            phase().boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
        }
        public void setDimensions(Space.Vector v) {dimensions.E(v); updateDimensions();}
        public double volume() {return dimensions.x * dimensions.y;}
        
        /** 
         * Computes origins for periodic images
         */
        public double[][] imageOrigins(int nShells) {
            int nImages = (2*nShells+1)*(2*nShells+1)-1;
            double[][] origins = new double[nImages][D];
            int k = 0;
            for(int i=-nShells; i<=nShells; i++) {
                for(int j=-nShells; j<=nShells; j++) {
                    if(i==0 && j==0) {continue;}
                    origins[k][0] = i*dimensions.x;
                    origins[k][1] = j*dimensions.y;
                    k++;
                }
            }
            return origins;
        }

        /** Returns coordinate shifts needed to draw all images that overflow into central image
         * 0, 1, or 3 shifts may be returned
         */
        int shiftX, shiftY;
        Vector r;
        public float[][] getOverflowShifts(Space.Vector rr, double distance) {
            shiftX = 0; shiftY = 0;
            r = (Vector)rr;
            
            if(r.x-distance < 0.0) {shiftX = +1;}
            else if(r.x+distance > dimensions.x) {shiftX = -1;}
            
            if(r.y-distance < 0.0) {shiftY = +1;}
            else if(r.y+distance > dimensions.y) {shiftY = -1;}
            
            if((shiftX == 0) && (shiftY == 0)) {
              shift = shift0;
            } else if((shiftX != 0) && (shiftY == 0)) {
              shift = new float[1][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
            } else if((shiftX == 0) && (shiftY != 0)) {
              shift = new float[1][D];
              shift[0][1] = (float)(shiftY*dimensions.y);
            } else if((shiftX != 0) && (shiftY != 0)) {
              shift = new float[3][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
              shift[1][1] = (float)(shiftY*dimensions.y);
              shift[2][0] = shift[0][0];
              shift[2][1] = shift[1][1];
            }
            return(shift);
        } //end of getOverflowShifts
    }  //end of BoundaryPeriodicSquare
    
    /**
     * Class for implementing slit periodic boundary conditions
     */
    protected static class BoundaryHSlit extends Boundary implements Space.Boundary.Periodic {
        public BoundaryHSlit() {this(Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryHSlit(Phase p) {this(p,Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryHSlit(Phase p, double lx, double ly) {super(p); dimensions.x = lx; dimensions.y = ly; updateDimensions();}
        public BoundaryHSlit(double lx, double ly) {dimensions.x = lx; dimensions.y = ly; updateDimensions();}
        public Space.Boundary.Type type() {return Boundary.HSLIT;}
        private final Vector temp = new Vector();
        private final Vector dimensions = new Vector();
        private final Vector dimensionsCopy = new Vector();
        private final Vector dimensionsHalf = new Vector();
        public final Space.Vector dimensions() {return dimensionsCopy;}
        private final void updateDimensions() {
            dimensionsHalf.Ea1Tv1(0.5,dimensions);
            dimensionsCopy.E(dimensions);
        }
        public Space.Vector randomPosition() {
            temp.x = dimensions.x*Simulation.random.nextDouble(); 
            temp.y = dimensions.y*Simulation.random.nextDouble(); 
            return temp;
        }
        public void nearestImage(Space.Vector dr) {nearestImage((Vector)dr);}
        public void nearestImage(Vector dr) {
           // dr.x -= dimensions.x * ((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x+0.5) : Math.ceil(dr.x/dimensions.x-0.5));
           // dr.y -= dimensions.y * ((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y+0.5) : Math.ceil(dr.y/dimensions.y-0.5));
           // final double dimxHalf = 0.5*dimensions.x;
           // final double dimyHalf = 0.5*dimensions.y;
            while(dr.y > +dimensionsHalf.y) {dr.y -= dimensions.y;}//System.out.println("Space2D here0 3");}
            while(dr.y < -dimensionsHalf.y) {dr.y += dimensions.y;}//System.out.println("Space2D here0 4");}       
        }
        public boolean centralImage(Coordinate c) {return centralImage(c.r);}
        public boolean centralImage(Space.Vector r) {return centralImage((Vector)r);}
        public boolean centralImage(Vector r) {
        /*    while(r.x > dimensions.x) r.x -= dimensions.x;
            while(r.x < 0.0)          r.x += dimensions.x;
            while(r.y > dimensions.y) r.y -= dimensions.y;
            while(r.y < 0.0)          r.y += dimensions.y;
        */  boolean changed = false;
        
            while(r.y > dimensions.y) {r.y -= dimensions.y; changed = true;}//System.out.println("Space2D here 3");}
            while(r.y < 0.0)          {r.y += dimensions.y; changed = true;}//System.out.println("Space2D here 4");}
            return changed;  
            //r.x -= dimensions.x * ((r.x >= 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
            //r.y -= dimensions.y * ((r.y >= 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y-1.0));
        }
        public void inflate(double scale) {
            dimensions.TE(scale); 
            updateDimensions();
            phase().boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
        }
        public void inflate(Space.Vector scale) {
            dimensions.TE(scale); 
            updateDimensions();
            phase().boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
        }
        public void setDimensions(Space.Vector v) {dimensions.E(v); updateDimensions();}
        public double volume() {return dimensions.x * dimensions.y;}
        
        /** 
         * Computes origins for periodic images
         */
        public double[][] imageOrigins(int nShells) {
            int nImages = (2*nShells+1)*(2*nShells+1)-1;
            double[][] origins = new double[nImages][D];
            int k = 0;
            for(int i=-nShells; i<=nShells; i++) {
                for(int j=-nShells; j<=nShells; j++) {
                    if(i==0 && j==0) {continue;}
                    origins[k][0] = i*dimensions.x;
                    origins[k][1] = j*dimensions.y;
                    k++;
                }
            }
            return origins;
        }

        /** Returns coordinate shifts needed to draw all images that overflow into central image
         * 0, 1, or 3 shifts may be returned
         */
        int shiftX, shiftY;
        Vector r;
        public float[][] getOverflowShifts(Space.Vector rr, double distance) {
            shiftX = 0; shiftY = 0;
            r = (Vector)rr;
            
            if(r.x-distance < 0.0) {shiftX = +1;}
            else if(r.x+distance > dimensions.x) {shiftX = -1;}
            
            if(r.y-distance < 0.0) {shiftY = +1;}
            else if(r.y+distance > dimensions.y) {shiftY = -1;}
            
            if((shiftX == 0) && (shiftY == 0)) {
              shift = shift0;
            } else if((shiftX != 0) && (shiftY == 0)) {
              shift = new float[1][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
            } else if((shiftX == 0) && (shiftY != 0)) {
              shift = new float[1][D];
              shift[0][1] = (float)(shiftY*dimensions.y);
            } else if((shiftX != 0) && (shiftY != 0)) {
              shift = new float[3][D];
              shift[0][0] = (float)(shiftX*dimensions.x);
              shift[1][1] = (float)(shiftY*dimensions.y);
              shift[2][0] = shift[0][0];
              shift[2][1] = shift[1][1];
            }
            return(shift);
        } //end of getOverflowShifts
    }  //end of BoundaryHSlit
  
  /*  public static final class BoundarySlidingBrick extends BoundaryPeriodicSquare {
        private double gamma = 0.0;
        private double delvx;
        private IntegratorMD.ChronoMeter timer;
        public BoundarySlidingBrick() {super();}
        public BoundarySlidingBrick(Phase p) {super(p);}
        public Space.Boundary.Type type() {return Boundary.SLIDING_BRICK;}
        public void setShearRate(double g) {gamma = g; computeDelvx();}
        public double getShearRate() {return gamma;}
        private void computeDelvx() {delvx = gamma*dimensions.y;}
        
        public void setTimer(IntegratorMD.ChronoMeter t) {timer = t;}
        
        public void nearestImage(Vector dr) {
            double delrx = delvx*timer.currentValue();
            double cory = ((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y+0.5) : Math.ceil(dr.y/dimensions.y-0.5));
            dr.x -= cory*delrx;
            dr.x -= dimensions.x * ((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x+0.5) : Math.ceil(dr.x/dimensions.x-0.5));
            dr.y -= dimensions.y * cory;
        }
        public void centralImage(Vector r) {
            double delrx = delvx*timer.currentValue();
            double cory = ((r.y >= 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y-1.0));
//            if(cory != 0.0) System.out.println(delrx*cory);
            r.x -= cory*delrx;
            r.x -= dimensions.x * ((r.x >= 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
            r.y -= dimensions.y * cory;
        }
        public void centralImage(Coordinate c) {
            Vector r = c.r;
            double cory = ((r.y > 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y-1.0));
            double corx = ((r.x > 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
            if(corx==0.0 && cory==0.0) return;
            double delrx = delvx*timer.currentValue();
            Vector p = c.p;
            r.x -= cory*delrx;
            r.x -= dimensions.x * corx; 
            r.y -= dimensions.y * cory;
            p.x -= cory*delvx;
        }
        
        public double[][] imageOrigins(int nShells) {
            int nImages = (2*nShells+1)*(2*nShells+1)-1;
            double[][] origins = new double[nImages][D];
            int k = 0;
            for(int i=-nShells; i<=nShells; i++) {
                for(int j=-nShells; j<=nShells; j++) {
                    if(i==0 && j==0) {continue;}
                    origins[k][0] = i*dimensions.x + j*delvx*timer.currentValue();
                    origins[k][1] = j*dimensions.y;
                    k++;
                }
            }
            return origins;
        }//end of imageOrigins
    }//end of BoundarySlidingBrick */
            
}//end of Space2D
