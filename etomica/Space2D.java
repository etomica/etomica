package etomica;
import java.awt.Graphics;
import java.awt.Color;
import etomica.units.*;

public class Space2D extends Space implements EtomicaElement {
    
    public static String version() {return "Space2D:01.07.07/"+Space.VERSION;}
    public static final int D = 2;
    public final int D() {return D;}
    
    public double sphereVolume(double r) {return Math.PI*r*r;}  //volume of a sphere of radius r
    public double sphereArea(double r) {return 2.0*Math.PI*r;}  //surface area of sphere of radius r (used for differential shell volume)
    public Space.Vector makeVector() {return new Vector();}
    public Space.Orientation makeOrientation() {return new Orientation();}
    public Space.Tensor makeTensor() {return new Tensor();}
    public Space.Coordinate makeCoordinate(Atom a) {//may want to revise this for o instanceof Molecule
        if(a instanceof AtomGroup) return new CoordinateGroup((AtomGroup)a);
        else if(a instanceof Atom && ((Atom)a).type instanceof AtomType.Rotator) return new OrientedCoordinate(a);
        else return new Coordinate(a);
    }
    public Space.CoordinatePair makeCoordinatePair(Phase p) {return new CoordinatePair(p);}
    
    public Space.Boundary.Type[] boundaryTypes() {return Boundary.TYPES;}
    public Space.Boundary makeBoundary() {return makeBoundary(Boundary.PERIODIC_SQUARE);}  //default
    public Space.Boundary makeBoundary(Space.Boundary.Type t) {
        if(t == Boundary.NONE) {return new BoundaryNone();}
        else if(t == Boundary.PERIODIC_SQUARE) {return new BoundaryPeriodicSquare();}
        else if(t == Boundary.SLIDING_BRICK) return new BoundarySlidingBrick();
        else return null;
    }
    
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
        public double x, y;
        public Vector () {x = 0.0; y = 0.0;}
        public Vector (double x, double y) {this.x = x; this.y = y;}
        public Vector (double[] a) {x = a[0]; y = a[1];}//should check length of a for exception
        public Vector (Vector u) {this.E(u);}
        public double[] toArray() {return new double[] {x, y};}
        public void sphericalCoordinates(double[] result) {
            result[0] = Math.sqrt(x*x + y*y);
            result[1] = Math.atan2(x, y);  //theta
        }
        public int length() {return D;}
        public int D() {return D;}
        public double component(int i) {return (i==0) ? x : y;}
        public void setComponent(int i, double d) {if(i==0) x=d; else y=d;}
        public void E(Vector u) {x = u.x; y = u.y;}
        public void E(double[] u) {x = u[0]; y = u[1];}  //should check length of array for exception
        public void E(double a) {x = a; y = a;}
//        public void E(int i, double a) {if(i==0) x = a; else y = a;}  //assumes i = 0 or 1
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
        public Space.Vector P(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x+u1.x; WORK.y = y+u1.y; return WORK;}
        public Space.Vector M(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x-u1.x; WORK.y = y-u1.y; return WORK;}
        public Space.Vector T(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x*u1.x; WORK.y = y*u1.y; return WORK;}
        public Space.Vector D(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x/u1.x; WORK.y = y/u1.y; return WORK;}
        public Space.Vector abs() {WORK.x = (x>0)?x:-x; WORK.y = (y>0)?y:-y; return WORK;}
        public double min() {return (x < y) ? x : y;}
        public double max() {return (x > y) ? x : y;}
        public double squared() {return x*x + y*y;}
        public double dot(Vector u) {return x*u.x + y*u.y;}
        public Space3D.Vector cross(Space3D.Vector u) {//not thread safe
            Space3D.Vector.WORK.x = y*u.z;
            Space3D.Vector.WORK.y = -x*u.z;
            Space3D.Vector.WORK.z = x*u.y - y*u.x;
            return Space3D.Vector.WORK;
        }
        public Space3D.Vector cross(Space2D.Vector u) {//not thread safe
            Space3D.Vector.WORK.x = 0.0;
            Space3D.Vector.WORK.y = 0.0;
            Space3D.Vector.WORK.z = x*u.y - y*u.x;
            return Space3D.Vector.WORK;
        }
        /**
         * Replaces this vector with its cross-product with the given 3D vector, with result projected
         * onto the 2D plane.  This vector becomes the result of (this vector) X u.
         */
        public void XE(Space3D.Vector u) {
            double xNew = y*u.z;
            y = -x*u.z;
            x = xNew;
        }
            
        public void normalize() {
            double norm = Math.sqrt(1/(x*x + y*y));
            x *= norm;
            y *= norm;
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
        public void randomDirection() {
            x = Math.cos(Math.PI*Simulation.random.nextDouble()); 
            y = Math.sqrt(1.0 - x*x);
            if(Simulation.random.nextDouble() < 0.5) y = -y;
        }
        public void E(Space.Vector u) {E((Vector)u);}
        public void PE(Space.Vector u) {PE((Vector)u);}
        public void ME(Space.Vector u) {ME((Vector)u);}
        public void TE(Space.Vector u) {TE((Vector)u);}
        public void DE(Space.Vector u) {DE((Vector)u);}
        public double dot(Space.Vector u) {return dot((Vector)u);}
    }
    
    public final static class Tensor extends Space.Tensor {
        public int length() {return D;}
        double xx, xy, yx, yy;
        public static final Tensor ZERO = new Tensor();
        public static final Tensor WORK = new Tensor();  //anything using WORK is not thread-safe
        public Tensor () {xx = xy = yx = yy = 0.0;}
        public Tensor (double xx, double xy, double yx, double yy) {this.xx = xx; this.xy = xy; this.yx = yx; this.yy = yy;}

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
    }

    public static final class CoordinatePair extends Space.CoordinatePair {  
        Coordinate c1;
        Coordinate c2;
        Boundary boundary;
        final Vector dimensions;   //assumes this is not transferred between phases
        private final Vector dr = new Vector();  //note that dr is not cloned if this is cloned -- should change this if using dr in clones; also this makes cloned coordinatePairs not thread-safe
        private double drx, dry, dvx, dvy;
        public CoordinatePair() {this(new BoundaryNone());}
        public CoordinatePair(Space.Boundary b) {boundary = (Boundary)b; dimensions = (Vector)boundary.dimensions();}
        public CoordinatePair(Phase p) {
            this(p.boundary());
            p.boundaryMonitor.addObserver(this);
        }
        /**
         * Implementation of Observer interface to update boundary if notified of change by phase.
         */
        public void update(java.util.Observable obs, Object arg) {boundary = (Boundary)arg;}
        public void reset(Space.Coordinate coord1, Space.Coordinate coord2) {  //don't usually use this; instead set c1 and c2 directly, without a cast
            c1 = (Coordinate)coord1;
            c2 = (Coordinate)coord2;
            reset();
        }
        public void reset() {
            dr.x = c2.r.x - c1.r.x;
            dr.y = c2.r.y - c1.r.y;
            boundary.nearestImage(dr);
            drx = dr.x; 
            dry = dr.y;
            r2 = drx*drx + dry*dry;
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
        public double position(int i) {return r.component(i);}
        public double momentum(int i) {return p.component(i);}
        public double kineticEnergy() {return 0.5*p.squared()*rm();}
        public void freeFlight(double t) {
            double tM = t*rm(); // t/mass
            r.x += p.x*tM;
            r.y += p.y*tM;
        }
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateBy(Space.Vector u) {r.PE((Space2D.Vector)u);}
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateBy(double d, Space.Vector u) {r.PEa1Tv1(d,(Vector)u);}
        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateTo(Space.Vector u) {r.E((Space2D.Vector)u);}      
        public void displaceBy(Space.Vector u) {rLast.E(r); translateBy((Vector)u);}
        public void displaceBy(double d, Space.Vector u) {rLast.E(r); translateBy(d,(Vector)u);}
        public void displaceTo(Space.Vector u) {rLast.E(r); translateTo((Vector)u);}  
        public void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        public void displaceToRandom(etomica.Phase p) {rLast.E(r); translateToRandom(p);}
        public void replace() {r.E(rLast);}
    //    public final void inflate(double s) {r.TE(s);}

        public void accelerateBy(Space.Vector u) {p.PE(u);}
        public void accelerateBy(double d, Space.Vector u) {p.PEa1Tv1(d,u);}

        public void randomizeMomentum(double temperature) {  //not very sophisticated; random only in direction, not magnitude
            double magnitude = Math.sqrt(mass()*temperature*(double)D);  //need to divide by sqrt(m) to get velocity
            momentum().setRandomSphere();
            momentum().TE(magnitude);
            //for debugging
      //      momentum().E(position());
      //      momentum().TE(magnitude/30.);
        }
    }
    
    public static class CoordinateGroup extends Coordinate implements Space.CoordinateGroup {
        public Coordinate firstChild, lastChild;
        public CoordinateGroup(AtomGroup a) {super(a);}

        public final Atom firstChild() {return (firstChild != null) ? firstChild.atom : null;}
        public final void setFirstChild(Atom a) {firstChild = (a != null) ? (Coordinate)a.coord : null;}
        public final Atom lastChild() {return (lastChild != null) ? lastChild.atom : null;}
        public final void setLastChild(Atom a) {lastChild = (a != null) ? (Coordinate)a.coord : null;}
        
        public Space.Vector position() {
            r.E(0.0); double massSum = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                r.PEa1Tv1(coord.mass(), coord.position()); massSum += coord.mass();
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
            double sum = 0.0;
            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                coord.freeFlight(t);
                if(coord == lastChild) break;
            }
        }
        public void translateBy(Space.Vector u) {
            Vector u0 = (Vector)u;
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
            work.E((Vector)u);
            work.ME(position());
            displaceBy(work);
        }
        public void displaceToRandom(etomica.Phase p) {
            displaceTo((Vector)p.boundary().randomPosition());
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
        public final void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        
        public void randomizeMomentum(double temperature) {
            switch(((AtomGroup)atom).childCount()) {
                case 0: return;
                case 1: firstChild.randomizeMomentum(temperature);//do not zero COM momentum if only one child atom
                        return;
                default://multi-atom group
                    work.E(0.0); double sum=0.0;
                    for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                        coord.randomizeMomentum(temperature);
                        work.PE(coord.momentum());
                        sum++;
                        if(coord == lastChild) break;
                    }
                    work.DE(-sum);
                    for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                        coord.accelerateBy(work);
                        if(coord == lastChild) break;
                    }
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
        public Space3D.Vector angularMomentum() {vector.z = L; return vector;}
        public Space3D.Vector angularVelocity() {vector.z = L/I[0]; return vector;}
        public void angularAccelerateBy(Space3D.Vector t) {L += t.z;}
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
        public static final Type SLIDING_BRICK = new Type("Sliding Brick");
        public static final Type[] TYPES = {NONE,PERIODIC_SQUARE,SLIDING_BRICK};
        public Boundary() {super();}
        public Boundary(Phase p) {super(p);}
        public abstract void nearestImage(Vector dr);
        public abstract void centralImage(Vector r);
        public abstract void centralImage(Coordinate c);
    }

    /**
     * Class for implementing no periodic boundary conditions
     */
    protected static class BoundaryNone extends Boundary {
        private final Vector temp = new Vector();
        private final double[][] shift0 = new double[0][D];
        public final Vector dimensions = new Vector();
        public Space.Vector dimensions() {return dimensions;}
        public BoundaryNone() {super();}
        public BoundaryNone(Phase p) {super(p);}
        public Space.Boundary.Type type() {return Boundary.NONE;}
        public final void nearestImage(Space.Vector dr) {}
        public final void centralImage(Space.Vector r) {}
        public final void nearestImage(Vector dr) {}
        public final void centralImage(Vector r) {}
        public final void centralImage(Coordinate c) {}
        public double volume() {return Double.MAX_VALUE;}
        public void inflate(double s) {}
        public double[][] imageOrigins(int nShells) {return new double[0][D];}
        public double[][] getOverflowShifts(Space.Vector rr, double distance) {return shift0;}
        public Space.Vector randomPosition() {  //arbitrary choice for this method in this boundary
            temp.x = Simulation.random.nextDouble(); 
            temp.y = Simulation.random.nextDouble(); 
            return temp;
        }
        public void draw(Graphics g, int[] origin, double scale) {}
    }//end of BoundaryNone

    /**
     * Class for implementing rectangular periodic boundary conditions
     */
    protected static class BoundaryPeriodicSquare extends Boundary implements Space.Boundary.Periodic {
        private final Vector temp = new Vector();
        private final double[][] shift0 = new double[0][D];
        private final double[][] shift1 = new double[1][D]; //used by getOverflowShifts
        private final double[][] shift3 = new double[3][D];
        public BoundaryPeriodicSquare() {this(Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p) {this(p,Default.BOX_SIZE,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p, double lx, double ly) {super(p); dimensions.x = lx; dimensions.y = ly;}
        public BoundaryPeriodicSquare(double lx, double ly) {dimensions.x = lx; dimensions.y = ly;}
        public Space.Boundary.Type type() {return Boundary.PERIODIC_SQUARE;}
        public final Vector dimensions = new Vector();
        public final Space.Vector dimensions() {return dimensions;}
        public Space.Vector randomPosition() {
            temp.x = dimensions.x*Simulation.random.nextDouble(); 
            temp.y = dimensions.y*Simulation.random.nextDouble(); 
            return temp;}
        public void nearestImage(Space.Vector dr) {nearestImage((Vector)dr);}
        public void nearestImage(Vector dr) {
            dr.x -= dimensions.x * ((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x+0.5) : Math.ceil(dr.x/dimensions.x-0.5));
            dr.y -= dimensions.y * ((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y+0.5) : Math.ceil(dr.y/dimensions.y-0.5));
        }
        public void centralImage(Coordinate c) {centralImage(c.r);}
        public void centralImage(Space.Vector r) {centralImage((Vector)r);}
        public void centralImage(Vector r) {
            r.x -= dimensions.x * ((r.x >= 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
            r.y -= dimensions.y * ((r.y >= 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y-1.0));
        }
        public void inflate(double scale) {dimensions.TE(scale);}
        public double volume() {return dimensions.x * dimensions.y;}
        
        public void draw(Graphics g, int[] origin, double scale) {
            g.setColor(Color.gray);
            double toPixels = scale*BaseUnit.Length.Sim.TO_PIXELS;
            g.drawRect(origin[0],origin[1],(int)(toPixels*dimensions.component(0))-1,(int)(toPixels*dimensions.component(1))-1);
            }
        /** Computes origins for periodic images
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
        public double[][] getOverflowShifts(Space.Vector rr, double distance) {
            Vector r = (Vector)rr;
            int shiftX = 0;
            int shiftY = 0;
            if(r.x-distance < 0.0) {shiftX = +1;}
            else if(r.x+distance > dimensions.x) {shiftX = -1;}
            
            if(r.y-distance < 0.0) {shiftY = +1;}
            else if(r.y+distance > dimensions.y) {shiftY = -1;}
            
            if(shiftX == 0) {
                if(shiftY == 0) {
                    return shift0;
                }
                else {
                    shift1[0][0] = 0.0;
                    shift1[0][1] = shiftY*dimensions.y;
                    return shift1;
                }
            }
            else { //shiftX != 0
                if(shiftY == 0) {
                    shift1[0][0] = shiftX*dimensions.x;
                    shift1[0][1] = 0.0;
                    return shift1;
                }
                else {
                    shift3[0][0] = shiftX*dimensions.x;
                    shift3[0][1] = 0.0;
                    shift3[1][0] = 0.0;
                    shift3[1][1] = shiftY*dimensions.y;
                    shift3[2][0] = shift3[0][0];
                    shift3[2][1] = shift3[1][1];
                    return shift3;
                }
            }
        } //end of getOverflowShifts
    }  //end of BoundaryPeriodicSquare
    
    public static final class BoundarySlidingBrick extends BoundaryPeriodicSquare {
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
    }//end of BoundarySlidingBrick
            
}//end of Space2D