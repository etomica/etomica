package etomica;
import java.awt.Graphics;
import java.awt.Color;
import java.util.Random;
import etomica.units.*;

public class Space1D extends Space implements EtomicaElement {
    
    public static final String version() {return "01.07.20";}
    
    public static final int D = 1;
    public static int drawingHeight = 10;  //height for drawing to 2D image
    public final int D() {return D;}
    
    public double sphereVolume(double r) {return 2.0*r;}  //volume of a sphere of radius r
    public double sphereArea(double r) {return 2.0;}      //surface area of sphere of radius r (used for differential shell volume)
    public Space.Vector makeVector() {return new Vector();}
    public Space.Orientation makeOrientation() {System.out.println("Orientation class not implemented in 1D"); return null;}
    public Space.Tensor makeTensor() {return new Tensor();}
    public Space.Coordinate makeCoordinate(Atom a) {return new Coordinate(a);}
    public Space.CoordinatePair makeCoordinatePair() {return new CoordinatePair();}

    public Space.Boundary.Type[] boundaryTypes() {return Boundary.TYPES;}
    public Space.Boundary makeBoundary() {return makeBoundary(Boundary.PERIODIC_SQUARE);}  //default
    public Space.Boundary makeBoundary(Space.Boundary.Type t) {
        if(t == Boundary.NONE) {return new BoundaryNone();}
        else if(t == Boundary.PERIODIC_SQUARE) {return new BoundaryPeriodicSquare();}
        else return null;
    }
    
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
        public static final Random random = new Random();
        public static final Vector ZERO = new Vector(0.0);
        public static final Vector WORK = new Vector();
        double x;
        public Vector () {x = 0.0;}
        public Vector (double a1) {x = a1;}
        public Vector (double[] a) {x = a[0];}//should check length of a for exception
        public Vector (Vector u) {this.E(u);}
        public String toString() {return "("+x+")";}
        public int length() {return D;}
        public int D() {return D;}
        public double component(int i) {return x;}
        public void setComponent(int i, double d) {x=d;}
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
        public Space.Vector P(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x+u1.x; return WORK;}
        public Space.Vector M(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x-u1.x; return WORK;}
        public Space.Vector T(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x*u1.x; return WORK;}
        public Space.Vector D(Space.Vector u) {Vector u1=(Vector)u; WORK.x = x/u1.x; return WORK;}
        public Space.Vector abs() {WORK.x = (x>0)?x:-x; return WORK;}
        public double min() {return x;}
        public double max() {return x;}
        public double squared() {return x*x;}
        public void normalize() {x = 1.0;}
        public double dot(Vector u) {return x*u.x;}
        public void randomStep(double d) {x += (2.*random.nextDouble()-1.0)*d;} //uniformly distributed random step in x and y, within +/- d
        public void setRandom(double d) {x = random.nextDouble()*d;}
        public void setRandom(double dx, double dy) {x = random.nextDouble()*dx;}
        public void setRandom(Vector u) {setRandom(u.x);}
        public void setRandomCube() {x = random.nextDouble() - 0.5;}
        public void setRandomSphere() {setRandomCube();}
        public void randomDirection() {x = (random.nextDouble() < 0.5) ? -1.0 : +1.0;}
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
    }
    
    public final static class Tensor extends Space.Tensor {
        double xx;
        public static final Tensor ZERO = new Tensor();
        public static final Tensor IDENTITY = new Tensor(1.0);
        public static final Tensor WORK = new Tensor();  //anything using WORK is not thread-safe
        public Tensor () {xx = 0.0;}
        public Tensor (double xx) {this.xx = xx;}

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
    }

    protected static final class CoordinatePair extends Space.CoordinatePair {  
        Coordinate c1;
        Coordinate c2;
        private final Vector dr = new Vector(); //note that dr is not cloned if this is cloned -- this should be changed if cloned vectors use dr; also this makes cloned coordinatePairs not thread-safe
        private double drx, dvx;
        public CoordinatePair() {super();}

        public void reset(Space.Coordinate coord1, Space.Coordinate coord2) {  //don't usually use this; instead set c1 and c2 directly, without a cast
            c1 = (Coordinate)coord1;
            c2 = (Coordinate)coord2;
            reset();
        }
        public void reset() {
            dr.x = c2.r.x - c1.r.x;
            c1.atom.parentPhase().boundary().nearestImage(dr);
            drx = dr.x; 
            r2 = drx*drx;
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

        public Space.Vector position() {return r;}
        public Space.Vector momentum() {return p;}
        public double position(int i) {return r.component(i);}
        public double momentum(int i) {return p.component(i);}
        public double kineticEnergy() {return 0.5*p.squared()*rm();}
        public void freeFlight(double t) {r.x += p.x*t*rm();}

        /**
        * Moves the atom by some vector distance
        * 
        * @param u
        */
        public void translateBy(Space.Vector u) {r.PE((Space1D.Vector)u);}
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
        public void translateTo(Space.Vector u) {r.E((Space1D.Vector)u);}      
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
        public abstract void centralImage(Vector r);
    }

    /**
     * Class for implementing no periodic boundary conditions
     */
    protected static final class BoundaryNone extends Boundary{
        private final Vector temp = new Vector();
        private final double[][] shift0 = new double[0][D];
        public final Vector dimensions = new Vector();
        public final Space.Vector dimensions() {return dimensions;}
        public static final Random random = new Random();
        public BoundaryNone() {super();}
        public BoundaryNone(Phase p) {super(p);}
        public Space.Boundary.Type type() {return Boundary.NONE;}
        public void nearestImage(Space.Vector dr) {}
        public void centralImage(Space.Vector r) {}
        public void nearestImage(Vector dr) {}
        public void centralImage(Vector r) {}
        public double volume() {return Double.MAX_VALUE;}
        public void inflate(double s) {}
        public double[][] imageOrigins(int nShells) {return new double[0][D];}
        public double[][] getOverflowShifts(Space.Vector rr, double distance) {return shift0;}
        public Space.Vector randomPosition() {  //arbitrary choice for this method in this boundary
            temp.x = random.nextDouble(); 
            return temp;
        }
        public void draw(Graphics g, int[] origin, double scale) {}
    }


    /**
     * Class for implementing rectangular periodic boundary conditions
     */
    protected static class BoundaryPeriodicSquare extends Boundary implements Space.Boundary.Periodic{
        private final Vector temp = new Vector();
        public static final Random random = new Random();
       //Explicit dimension to 2 because drawing to 1D image
        private final double[][] shift0 = new double[0][2];
        private final double[][] shift1 = new double[1][2]; //used by getOverflowShifts
        public BoundaryPeriodicSquare() {this(Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p) {this(p,Default.BOX_SIZE);}
        public BoundaryPeriodicSquare(Phase p, double lx) {super(p); dimensions.x = lx;}
        public BoundaryPeriodicSquare(double lx) {dimensions.x = lx;}
        public Space.Boundary.Type type() {return Boundary.PERIODIC_SQUARE;}
        public final Vector dimensions = new Vector();
        public final Space.Vector dimensions() {return dimensions;}
        public Space.Vector randomPosition() {
            temp.x = dimensions.x*random.nextDouble(); 
            return temp;}
        public void nearestImage(Space.Vector dr) {nearestImage((Vector)dr);}
        public void nearestImage(Vector dr) {
            dr.x -= dimensions.x * ((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x+0.5) : Math.ceil(dr.x/dimensions.x-0.5));
        }
        public void centralImage(Space.Vector r) {centralImage((Vector)r);}
        public void centralImage(Vector r) {
            r.x -= dimensions.x * ((r.x >= 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
        }
        public void inflate(double scale) {dimensions.TE(scale);}
        public double volume() {return dimensions.x;}
        public void draw(Graphics g, int[] origin, double scale) {
            g.setColor(Color.gray);
            double toPixels = scale*BaseUnit.Length.Sim.TO_PIXELS;
            g.drawRect(origin[0],origin[1],(int)(toPixels*dimensions.component(0))-1,(int)(toPixels*dimensions.component(1))-1);
            }
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
        public double[][] getOverflowShifts(Space.Vector rr, double distance) {
            Vector r = (Vector)rr;
            int shiftX = 0;
            if(r.x-distance < 0.0) {shiftX = +1;}
            else if(r.x+distance > dimensions.x) {shiftX = -1;}
            
            if(shiftX == 0) {
                return shift0;
            }
            else { //shiftX != 0
                shift1[0][0] = shiftX*dimensions.x;
                return shift1;
            }
        } //end of getOverflowShifts
    }  //end of BoundaryPeriodic
            
}//end of Space1D