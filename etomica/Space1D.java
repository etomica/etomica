package etomica;

import etomica.statmech.MaxwellBoltzmann;

//CoordinateGroup is not updated to the same structure as used in Space2D and Space3D
//centralImage not updated to molecule form as in Space3D

 /* History of changes
  * 09/01/02 (DAK) added accelerateTo method to Coordinate
  *                changed CoordinateGroup.randomizeMomentum to not enforce zero COM momentum
  * 09/05/02 (DAK) fixed error in accelerateTo (still probably does not do what one expects
  *                if accelerating to nonzero momentum).
  * 07/10/03 (DAK) added resetV method to CoordinatePair
  * 08/27/03 (DAK) added isZero method to Vector
  * 08/29/03 (DAK) implemented centralImage(Space.Coordinate) in Boundary
  * 12/09/03 (DAK) changed setRandomSphere to give point on surface of sphere
  * (just at +/- 1 instead of between +/- 1); added setRandomInSphere method
  * 01/22/04 (DAK) added positionCOM and translateCOMTo to CoordinateGroup;
  * redefined position() in CoordinateGroup to be first-atom position (as it has
  * been for Space3D for some time now).
  */
public class Space1D extends Space implements EtomicaElement {
    
    public static final int D = 1;
    public static int drawingHeight = 10;  //height for drawing to 2D image
    public final int D() {return D;}
    public final int powerD(int n) {return n;}
    public final double powerD(double a) {return a;}
    public static final Vector ORIGIN = new Vector();
    public final Space.Vector origin() {return ORIGIN;}
    public static final Space1D INSTANCE = new Space1D();
    
    public Space1D() {super(1);}
    
    public double sphereVolume(double r) {return 2.0*r;}  //volume of a sphere of radius r
    public double sphereArea(double r) {return 2.0;}      //surface area of sphere of radius r (used for differential shell volume)
    public Space.Vector makeVector() {return new Vector();}
    public Space.Orientation makeOrientation() {System.out.println("Orientation class not implemented in 1D"); return null;}
    public Space.Tensor makeTensor() {return new Tensor();}
    public Space.Tensor makeRotationTensor() {return new RotationTensor();}
    public Space.Coordinate makeCoordinate(Atom a) {
        if(a.node instanceof AtomTreeNodeGroup) return new CoordinateGroup(a);
//        else if(a.type instanceof AtomType.Rotator) return new OrientedCoordinate(a);
        return new Coordinate(a);
    }
    public Space.CoordinatePair makeCoordinatePair() {return new CoordinatePair();}

    public Space.Boundary.Type[] boundaryTypes() {return Boundary.TYPES;}
    public Space.Boundary makeBoundary() {return makeBoundary(Boundary.PERIODIC_SQUARE);}  //default
    public Space.Boundary makeBoundary(Space.Boundary.Type t) {
        if(t == Boundary.NONE) {return new BoundaryNone();}
        else if(t == Boundary.PERIODIC_SQUARE) {return new BoundaryPeriodicSquare();}
        else return null;
    }
    
    public int[] makeArrayD(int i) {return new int[] {i};}
    public double[] makeArrayD(double d) {return new double[] {d};}
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("One-dimensional space");
        return info;
    }

    public static final double r2(Vector u1, Vector u2, Boundary b) {
        Vector.WORK.x = u1.x - u2.x;
        b.nearestImage(Vector.WORK);
        return Vector.WORK.x*Vector.WORK.x;
    }
        
    public final static class Vector extends Space.Vector {  //declared final for efficient method calls
        public static final Vector ZERO = new Vector(0.0);
        public static final Vector WORK = new Vector();
        private double x;
        public Vector () {x = 0.0;}
        public Vector (double a1) {x = a1;}
        public Vector (double[] a) {x = a[0];}//should check length of a for exception
        public Vector (Vector u) {this.E(u);}
        public boolean equals(Space.Vector v) {return equals((Vector)v);}
        public boolean equals(Vector v) {return (x == v.x);}
        public boolean isZero() {return x==0;}
        public String toString() {return "("+x+")";}
        public int length() {return D;}
        public int D() {return D;}
        public double x(int i) {return x;}
        public void setX(int i, double d) {x=d;}
        public double[] toArray() {return new double[] {x};}
        public void sphericalCoordinates(double[] result) {result[0] = x;}
        public void E(Vector u) {x = u.x;}
        public void E(double a) {x = a;}
        public void E(double[] a) {x = a[0];}
        public void E(int i, double a) {x = a;}  //assumes i = 0
        public void Ea1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x = a1*u1.x;}
        public void PEa1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x += a1*u1.x;}
        public void PE(Vector u) {x += u.x;}
        public void PE(double a) {x += a;}
        public void ME(Vector u) {x -= u.x;}
        public void PE(int i, double a) {x += a;}
        public void TE(double a) {x *= a;}
        public void TE(Vector u) {x *= u.x;}
        public void TE(int i, double a) {x *= a;}
        public void DE(double a) {x /= a;}
        public void DE(Vector u) {x /= u.x;}
        public void Ev1Pv2(Space.Vector u1, Space.Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x + v2.x;
        }
        public void Ev1Mv2(Space.Vector u1, Space.Vector u2) {
            Vector v1 = (Vector)u1; Vector v2 = (Vector)u2;
            x = v1.x - v2.x;
        }
		public double Mv1Squared(Space.Vector u1) {
			double dx = x-((Vector)u1).x;
			return dx*dx;
		}
        public void mod(Space.Vector u) {
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


        public Space.Vector P(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x+u1.x; return WORK;}
        public Space.Vector M(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x-u1.x; return WORK;}
        public Space.Vector T(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x*u1.x; return WORK;}
        public Space.Vector D(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x/u1.x; return WORK;}
        public void abs() {x = (x>0)?x:-x;}
        public double min() {return x;}
        public double max() {return x;}
        public double squared() {return x*x;}
        public void normalize() {x = 1.0;}
        public double dot(Vector u) {return x*u.x;}
        public void transform(Space.Tensor A) {transform((Tensor)A);}
        public void transform(Tensor A) {x = A.xx * x;}
        public void transform(Space.Boundary b, Space.Vector r0, Space.Tensor A) {transform((Boundary)b, (Vector)r0, (Tensor)A);}
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
        public void E(Space.Vector u) {E((Vector)u);}
        public void PE(Space.Vector u) {PE((Vector)u);}
        public void ME(Space.Vector u) {ME((Vector)u);}
        public void TE(Space.Vector u) {TE((Vector)u);}
        public void DE(Space.Vector u) {DE((Vector)u);}
        public double dot(Space.Vector u) {return dot((Vector)u);}
        /**
         * Sets this vector equal to its cross product of with a 3D vector.
         * Result is projected into this space, and thus outcome is to make this vector zero.
         */
        public void XE(Space3D.Vector u) {x = 0.0;}
        public Space3D.Vector cross(Space2D.Vector u) {return null;}
        public Space3D.Vector cross(Space3D.Vector u) {return null;}
        public void randomRotate(double deltheta){//no implementation in 1D
        }
    }
    
    public static class Tensor implements Space.Tensor {
        double xx;
        public static final Tensor ZERO = new Tensor();
        public static final Tensor IDENTITY = new Tensor(new double[] {1.0});
        public static final Tensor WORK = new Tensor();  //anything using WORK is not thread-safe
        public Tensor () {xx = 0.0;}
        public Tensor (double[] d) {this.E(d);}

        public int length() {return D;}
        public double component(int i, int j) {return xx;}
        public void setComponent(int i, int j, double d) {xx = d;}
        public void E(Tensor t) {xx=t.xx;}
        public void E(Vector u1, Vector u2) {xx=u1.x*u2.x;}
        public void E(double a) {xx = a;}
        public void PE(Tensor t) {xx+=t.xx;}
        public void PE(int i, int j, double a) {xx += a;}
        public void PE(Vector u1, Vector u2) {xx+=u1.x*u2.x;}
        public double trace() {return xx;}
        
        public void E(Space.Tensor t) {E((Tensor)t);}
        public void E(Space.Vector u1, Space.Vector u2) {E((Vector)u1, (Vector)u2);}
        public void PE(Space.Tensor t) {PE((Tensor)t);}
        public void PE(Space.Vector u1, Space.Vector u2) {PE((Vector)u1, (Vector)u2);}
        public void TE(double a) {xx*=a;}
        public void E(double[] d) {
            if(d.length != 1) throw new IllegalArgumentException("Array size incorrector for tensor");
            xx = d[0];
        }
        public void assignTo(double[] d) {
            if(d.length != 1) throw new IllegalArgumentException("Array size incorrector for tensor");
            d[0] = xx;
        }
        
    }

    public static class RotationTensor extends Tensor implements Space.RotationTensor {
        public RotationTensor() {super(); reset();}
        public void reset() {
            xx = 1.0;
        }
        public void setAxial(int i, double theta) {
        }
        public void setAngles(double[] angles) {}
        public void invert() {}
    }

    protected static final class CoordinatePair extends Space.CoordinatePair {
        Coordinate c1;
        Coordinate c2;
        private final Vector dr = new Vector(); //note that dr is not cloned if this is cloned -- this should be changed if cloned vectors use dr; also this makes cloned coordinatePairs not thread-safe
        private double drx, dvx;
		private NearestImageTransformer nearestImageTransformer = Space.Boundary.NULL;

		public void setNearestImageTransformer(NearestImageTransformer b) {this.nearestImageTransformer = b;}
		public NearestImageTransformer getNearestImageTransformer() {return nearestImageTransformer;}		
        public double r2() {return drx*drx;}
        public void reset(Space.Coordinate coord1, Space.Coordinate coord2) {  //don't usually use this; instead set c1 and c2 directly, without a cast
            c1 = (Coordinate)coord1;
            c2 = (Coordinate)coord2;
            reset();
        }
        public void trueReset(Space.Coordinate coord1, Space.Coordinate coord2, double falseTime) {
            c1 = (Coordinate)coord1;
            c2 = (Coordinate)coord2;
            trueReset(falseTime);
        }
        public void reset() {
            dr.x = c2.r.x - c1.r.x;
            c1.atom.node.parentPhase().boundary().nearestImage(dr);
            drx = dr.x;
        }
        public void trueReset(double falseTime) {
            resetV();
            dr.Ev1Mv2(c2.r,c1.r);

            dr.x += falseTime * dvx;
            nearestImageTransformer.nearestImage(dr);
            drx = dr.x;
        }
        public void resetV() {
            double rm1 = c1.rm();
            double rm2 = c2.rm();
            dvx = (rm2*c2.p.x - rm1*c1.p.x);  
        }
        /**
         * Recomputes pair separation, with atom 2 shifted by the given vector
         * Does not apply any PBC, regardless of boundary chosen for space
         */
        public void reset(Space1D.Vector M) {
            dr.x = c2.r.x - c1.r.x + M.x;
            drx = dr.x;
            r2 = drx*drx;
        }
        public Space.Vector dr() {return dr;}
        public double dr(int i) {return drx;}
        public double dv(int i) {return dvx;}
        public double v2() {
            return dvx*dvx;
        }
        public double vDot(Space.Vector u) {return vDot((Space1D.Vector)u);}
        public double vDot(Space1D.Vector u) {return dvx*u.x;}
        public double vDotr() {
            return drx*dvx;
        }
        public void push(double impulse) {  //changes momentum in the direction joining the atoms
            c1.p.x += impulse*drx;
            c2.p.x -= impulse*drx;
        }
        public void truePush(Space.Vector u, double falseTime) {
            c1.p.PE(u);
            c2.p.ME(u);
            
            c1.r.PEa1Tv1(-falseTime,u);
            c2.r.PEa1Tv1(falseTime,u);
        }
        public void nudge(double rDelta) {
            double ratio = c2.mass()*c1.rm()*rDelta;
            c1.r.x -= ratio*dr.x;
            c2.r.x += ratio*dr.x;
        }
        public void setSeparation(double r2New) {
            double ratio = c2.mass()*c1.rm();  // (mass2/mass1)
            double delta = (Math.sqrt(r2New/this.r2()) - 1.0)/(1+ratio);
            c1.r.x -= ratio*delta*drx;
            c2.r.x += delta*drx;
            //need call reset?
        }
    }

    public static class Coordinate extends Space.Coordinate {
        public Coordinate nextCoordinate, previousCoordinate;
        public final Vector r = new Vector();  //Cartesian coordinates
        public final Vector p = new Vector();  //Momentum vector
        public final Vector rLast = new Vector();  //vector for saving position
        protected final Vector work = new Vector();
        public Coordinate(Atom atom) {super(atom);}
        
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

        public void transform(Space.Vector r0, Space.Tensor A) {
            r.transform((Boundary)atom.node.parentPhase().boundary(), (Vector)r0, (Tensor)A);
        }
        public Space.Vector position() {return r;}
        public Space.Vector truePosition(double falseTime) {
            work.E(r);
            work.PEa1Tv1(falseTime*rm(),p);
            return work;
        }
        public Space.Vector momentum() {return p;}
        public double position(int i) {return r.x(i);}
        public double truePosition(int i, double falseTime) {return r.x(i)+falseTime*rm()*p.x(i);}
        public double momentum(int i) {return p.x(i);}
        public double kineticEnergy() {return 0.5*p.squared()*rm();}
        public void freeFlight(double t) {r.x += p.x*t*rm();}
        public void inflate(double s) {r.x *= s;}
        public void inflate(Space.Vector s) {r.x *= ((Vector)s).x;}

        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateBy(Space.Vector u) {
            r.PE((Space1D.Vector)u);
        }
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateBy(double d, Space.Vector u) {
            r.PEa1Tv1(d,u);
        }
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateTo(Space.Vector u) {
            r.E((Space1D.Vector)u);
        }      
        public void replace() {
            r.E(rLast);
        }
        public void displaceBy(Space.Vector u) {rLast.E(r); translateBy(u);}
        public void displaceBy(double d, Space.Vector u) {rLast.E(r); translateBy(d,u);}
        public void displaceTo(Space.Vector u) {rLast.E(r); translateTo(u);}  
        public void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        public void displaceToRandom(etomica.Phase phase) {rLast.E(r); translateToRandom(phase);}
    //    public final void inflate(double s) {r.TE(s);}

        public void accelerateBy(Space.Vector u) {p.PE(u);}
        public void accelerateBy(double d, Space.Vector u) {p.PEa1Tv1(d,u);}
        public void accelerateTo(Space.Vector u) {p.E(u);}
        public void trueAccelerateTo(Space.Vector u, double falseTime) {
            r.x -= falseTime * (((Vector)u).x - p.x);
            p.x = ((Vector)u).x;
        }
        
        public void randomizeMomentum(double temperature) {
            if(isStationary()) {p.E(0.0); return;}
            p.setX(0,MaxwellBoltzmann.randomMomentumComponent(temperature,mass()));
        }
    }
    
    public static class CoordinateGroup extends Coordinate {
        public Coordinate firstChild, lastChild;//to be removed
        public Atom firstAtom;
        public CoordinateGroup(Atom a) {super(a);}

		//these may be removed at some point, to use childIterator instead (as in other space classes)
        public final Atom firstAtom() {return (firstChild != null) ? firstChild.atom : null;}
        public final void setFirstAtom(Atom a) {firstChild = (a != null) ? (Coordinate)a.coord : null;}
        public final Atom lastAtom() {return (lastChild != null) ? lastChild.atom : null;}
        public final void setLastAtom(Atom a) {lastChild = (a != null) ? (Coordinate)a.coord : null;}
        
        /**
         * Applies transformation to COM of group, keeping all internal atoms at same relative
         * positions.
         */
        public void transform(Space.Vector r0, Space.Tensor A) {
            work.E(position()); //work = r
            work.transform((Boundary)atom.node.parentPhase().boundary(), (Vector)r0, (Tensor)A);
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
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                r.PEa1Tv1(coord.mass(), coord.positionCOM()); massSum += coord.mass();
                if(coord == lastChild) break;
            }
            r.DE(massSum);
            return r;
        }
        public Space.Vector momentum() {
            p.E(0.0);
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                p.PE(coord.momentum());
                if(coord == lastChild) break;
            }
            return p;
        }
        public double position(int i) {
            double sum = 0.0; double massSum = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                sum += coord.mass()*coord.position(i); massSum += coord.mass();
                if(coord == lastChild) break;
            }
            sum /= massSum;
            return sum;
        }
        public double momentum(int i) {
            double sum = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                sum += coord.mass()*coord.momentum(i);
                if(coord == lastChild) break;
            }
            return sum;
        }
        public double kineticEnergy() {
            double sum = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                sum += coord.kineticEnergy();
                if(coord == lastChild) break;
            }
            return sum;
        }
        public void freeFlight(double t) {
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.freeFlight(t);
                if(coord == lastChild) break;
            }
        }
        public void inflate(double scale) {
            work.E(position());
            work.TE(scale-1.0);
            displaceBy(work);
        }
        public void inflate(Space.Vector scale) {
            scale.PE(-1.0);
            work.E(position());
            work.TE(scale);
            displaceBy(work);
            scale.PE(1.0);
        }
        public void translateBy(Space.Vector u) {translateBy((Vector)u);}
        public void translateBy(Vector u0) {
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.translateBy(u0);
                if(coord == lastChild) break;
            }
        }
        public void translateBy(double d, Space.Vector u) {
            Vector u0 = (Vector)u;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.translateBy(d, u0);
                if(coord == lastChild) break;
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
            work.Ea1Tv1(-1,position()); //position() uses work, so need this first
            work.PE((Vector)u);
            translateBy(work);
        }
        public void displaceBy(Space.Vector u) {
            Vector u0 = (Vector)u;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.displaceBy(u0);
                if(coord == lastChild) break;
            }
        }
        public void displaceBy(double d, Space.Vector u) {
            Vector u0 = (Vector)u;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.displaceBy(d, u0);
                if(coord == lastChild) break;
            }
        }
        public void displaceTo(Space.Vector u) {
            work.Ea1Tv1(-1,position()); //position() uses work, so need this first
            work.PE((Vector)u);
            displaceBy(work);
        }
        public void displaceToRandom(etomica.Phase phase) {
            displaceTo(phase.boundary().randomPosition());
        }
        public void replace() {
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.replace();
                if(coord == lastChild) break;
            }
        }
        public void accelerateBy(Space.Vector u) {
            Vector u0 = (Vector)u;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.accelerateBy(u0);
                if(coord == lastChild) break;
            }
        }
        public void accelerateBy(double d, Space.Vector u) {
            Vector u0 = (Vector)u;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.accelerateBy(d, u0);
                if(coord == lastChild) break;
            }
        }
        public void accelerateTo(Space.Vector u) {
    //        throw new RuntimeException("Space1D.CoordinateGroup.accelerateTo not implemented");
            //ugly
            int sum = 0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                sum ++;
                if(coord == lastChild) break;
            }

            work.Ea1Tv1(-1.0/(double)sum,momentum());//probably need this first
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
        //            work.E(0.0); double sum=0.0;
                    for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                        coord.randomizeMomentum(temperature);
        //                work.PE(coord.momentum());
        //                sum++;
        //                if(coord == lastChild) break;
                    }
        //            work.DE(-sum);
        //            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
        //                coord.accelerateBy(work);
        //                if(coord == lastChild) break;
        //            }
            }//end switch
        }//end randomizeMomentum
    }//end of CoordinateGroup
    
    public static abstract class Boundary extends Space.Boundary {
        public static class Type extends Space.Boundary.Type {
            private Type(String label) {super(label);}
            public Constants.TypedConstant[] choices() {return TYPES;}
        }
        public static final Type NONE = new Type("None");
        public static final Type PERIODIC_SQUARE = new Type("Periodic");
        public static final Type HARMONIC = new Type("Harmonic");
        public static final Type HARD = new Type("Hard");
        public static final Type[] TYPES = {NONE, PERIODIC_SQUARE ,HARMONIC, HARD};
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
    protected static final class BoundaryNone extends Boundary{
        private final Vector temp = new Vector();
        private final Vector dimensions = new Vector(Default.BOX_SIZE);
        private final Vector dimensionsCopy = new Vector();
        public final Space.Vector dimensions() {dimensionsCopy.E(dimensions); return dimensionsCopy;}
        public BoundaryNone() {super();}
        public BoundaryNone(Phase p) {super(p);}
        public Space.Boundary.Type type() {return Boundary.NONE;}
        public void nearestImage(Space.Vector dr) {}
		public boolean centralImage(Space.Vector r) {return false;}
		public boolean centralImage(Coordinate c) {return false;}
        public void nearestImage(Vector dr) {}
        public boolean centralImage(Vector r) {return false;}
        public double volume() {return dimensions.x;}
        public void inflate(double s) {dimensions.TE(s);}
        public void inflate(Space.Vector s) {dimensions.TE(s);}
        public void setDimensions(Space.Vector v) {dimensions.E(v);}
        public double[][] imageOrigins(int nShells) {return new double[0][D];}
        public float[][] getOverflowShifts(Space.Vector rr, double distance) {return shift0;}
        public Space.Vector randomPosition() {  //arbitrary choice for this method in this boundary
            temp.x = dimensions.x*Simulation.random.nextDouble(); 
            return temp;
        }
    }


    /**
     * Class for implementing rectangular periodic boundary conditions
     */
    protected static class BoundaryPeriodicSquare extends Boundary implements Space.Boundary.Periodic{
        public BoundaryPeriodicSquare() {this(Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p) {this(p,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p, double lx) {super(p); dimensions.x = lx; updateDimensions();}
        public BoundaryPeriodicSquare(double lx) {dimensions.x = lx; updateDimensions();}
        public Space.Boundary.Type type() {return Boundary.PERIODIC_SQUARE;}
        private final Vector temp = new Vector();
        private final Vector dimensions = new Vector();
        private final Vector dimensionsHalf = new Vector();
        private final Vector dimensionsCopy = new Vector();
        private final void updateDimensions() {
            dimensionsHalf.Ea1Tv1(0.5,dimensions);
            dimensionsCopy.E(dimensions);
        }
        public final Space.Vector dimensions() {return dimensionsCopy;}
        public Space.Vector randomPosition() {
            temp.x = dimensions.x*Simulation.random.nextDouble(); 
            return temp;
        }
        public void nearestImage(Space.Vector dr) {nearestImage((Vector)dr);}
        public void nearestImage(Vector dr) {
            while(dr.x > +dimensionsHalf.x) dr.x -= dimensions.x;
            while(dr.x < -dimensionsHalf.x) dr.x += dimensions.x;
        }
        public boolean centralImage(Coordinate c) {return centralImage(c.position());}
        public boolean centralImage(Space.Vector r) {return centralImage((Vector)r);}
        public boolean centralImage(Vector r) {
            boolean changed = false;
            while(r.x > dimensions.x) {r.x -= dimensions.x; changed = true;}
            while(r.x < 0.0)          {r.x += dimensions.x; changed = true;}
            return changed;
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
        public double volume() {return dimensions.x;}
        /** Computes origins for periodic images
         */
        public double[][] imageOrigins(int nShells) {
            int nImages = 2*nShells;
            double[][] origins = new double[nImages][D];
            int k = 0;
            for(int i=-nShells; i<=nShells; i++) {
                if(i==0) {continue;}
                origins[k][0] = i*dimensions.x;
                k++;
            }
            return origins;
        }

        /** Returns coordinate shifts needed to draw all images that overflow into central image
         * 0, or 1 shifts may be returned
         */
        int shiftX;
        Vector r;
        public float[][] getOverflowShifts(Space.Vector rr, double distance) {
            r = (Vector)rr;
            shiftX = 0;
            if(r.x-distance < 0.0) {shiftX = +1;}
            else if(r.x+distance > dimensions.x) {shiftX = -1;}
            
            if(shiftX == 0)
                return shift0;
            shift = new float[1][2];
            shift[0][0] = (float)(shiftX*dimensions.x);
            return shift;
        } //end of getOverflowShifts
    }  //end of BoundaryPeriodic
    
    
//    public static void main(String[] args) {
//        Default.ATOM_SIZE = 1.0;
//        etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic(new Space1D());
//  //    setIteratorFactory(new IteratorFactoryCell(this));
//        AtomIteratorNeighbor nbrIterator = new AtomIteratorNeighbor();
//	    IntegratorMC integrator = new IntegratorMC(sim);
//	    MCMoveAtom mcmoveAtom = new MCMoveAtom(integrator);
//	    Species species = new SpeciesSpheresMono(sim);
//	    
//	    species.setNMolecules(3);
//	    final Phase phase = new Phase(sim);
//	    Potential2 potential = new P2HardSphere();
////	    Potential2 potentialOrder = new P2XOrder();
////	    potential.setSpecies(species,species);
////	    potentialOrder.setSpecies(species,species);
//	    Controller controller = new Controller(sim);
//	    etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase(sim);
//        display.setColorScheme(new etomica.graphics.ColorSchemeRandom());
//		sim.elementCoordinator.go();
//		
//		AtomPairIterator api = new ApiGeneral(sim.space, new AtomIteratorList(), nbrIterator);
//		potential.setIterator(api);
//		
//		NeighborManager.Criterion criterion = new NeighborManager.Criterion() {
//		    public boolean areNeighbors(Atom a1, Atom a2) {
//		        return Math.abs(a1.node.index()-a2.node.index()) == 1 ||
//		           (a1==phase.firstAtom() && a2==phase.lastAtom());
//		    }};
//		nbrIterator.setupNeighbors(phase.speciesMaster.atomList, criterion);
//        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
//    }
            
}//end of Space1D