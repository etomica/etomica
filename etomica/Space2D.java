package simulate;
import java.awt.Graphics;
import java.util.Random;

public class Space2D extends Space {
    
    public static final int D = 2;
    public final int D() {return D;}

    public Space.AtomCoordinate makeAtomCoordinate(Atom a) {return new AtomCoordinate(a);}
    public Space.Coordinate makeCoordinate() {return new Coordinate();}
    public Space.Vector makeVector() {return new Vector();}
    public final simulate.AtomPair.Iterator.A makePairIteratorFull(Space.Boundary boundary, Atom iF, Atom iL, Atom oF, Atom oL) {return new PairIteratorFull((Boundary)boundary,iF,iL,oF,oL);}
    public final simulate.AtomPair.Iterator.A makePairIteratorHalf(Space.Boundary boundary, Atom iL, Atom oF, Atom oL) {return new PairIteratorHalf((Boundary)boundary,iL,oF,oL);}
    public final simulate.AtomPair.Iterator.A makePairIteratorFull(Space.Boundary boundary) {return new PairIteratorFull((Boundary)boundary);}
    public final simulate.AtomPair.Iterator.A makePairIteratorHalf(Space.Boundary boundary) {return new PairIteratorHalf((Boundary)boundary);}
    public simulate.AtomPair makeAtomPair(Space.Boundary boundary, Atom a1, Atom a2) {return new AtomPair((Boundary)boundary, a1, a2);}
    public Space.Boundary makeBoundary(int b) {
        switch(b) {
            case(Boundary.NONE):            return new BoundaryNone();
            case(Boundary.PERIODIC_SQUARE): return new BoundaryPeriodicSquare();
            default:                        return null;
        }
    }
        
    public static abstract class Boundary implements Space.Boundary {
        public static final int NONE = 0;
        public static final int PERIODIC_SQUARE = 1;
        public final Vector dimensions = new Vector();
        public void centralImage(Space.Vector r) {centralImage((Vector)r);}
        public final Space.Vector dimensions() {return dimensions;}
        public abstract void centralImage(Vector r);
        public abstract void apply(Vector r);
        public abstract double volume();
        public abstract void inflate(double s);
    }

    /**
     * Class for implementing no periodic boundary conditions
     */
    private static final class BoundaryNone extends Boundary {
        private final Vector temp = new Vector();
        private final double[][] shift0 = new double[0][D];
        public static final Random random = new Random();
        public void apply(Vector dr) {}
        public void centralImage(Vector r) {}
        public double volume() {return Double.MAX_VALUE;}
        public void inflate(double s) {}
        public double[][] imageOrigins(int nShells) {return new double[0][D];}
        public double[][] getOverflowShifts(Space.Vector rr, double distance) {return shift0;}
        public Space.Vector randomPosition() {temp.x = random.nextDouble()-0.5; temp.y = random.nextDouble()-0.5; return temp;}
    }

    /**
     * Class for implementing rectangular periodic boundary conditions
     */
    private static final class BoundaryPeriodicSquare extends Boundary {
        private final Vector temp = new Vector();
        public static final Random random = new Random();
        private final double[][] shift0 = new double[0][D];
        private final double[][] shift1 = new double[1][D]; //used by getOverflowShifts
        private final double[][] shift3 = new double[3][D];
        public BoundaryPeriodicSquare() {this(1.0,1.0);}
        public BoundaryPeriodicSquare(double lx, double ly) {dimensions.x = lx; dimensions.y = ly;}
        public Space.Vector randomPosition() {temp.x = random.nextDouble()-0.5; temp.y = random.nextDouble()-0.5; return temp;}
        public void apply(Vector dr) {
            dr.x -= (dr.x > 0.0) ? Math.floor(dr.x+0.5) : Math.ceil(dr.x-0.5);
            dr.y -= (dr.y > 0.0) ? Math.floor(dr.y+0.5) : Math.ceil(dr.y-0.5);
//            if(dr.x > 0) {dr.x -= Math.floor(dr.x+0.5);}
//            else         {dr.x -= Math.ceil(dr.x-0.5);}
//            if(dr.y > 0) {dr.y -= Math.floor(dr.y+0.5);}
//            else         {dr.y -= Math.ceil(dr.y-0.5);}
            dr.x *= dimensions.x;
            dr.y *= dimensions.y;
        }
        public void centralImage(Space.Vector r) {centralImage((Vector)r);}
        public void centralImage(Vector r) {
            r.x -= (r.x > 0.0) ? Math.floor(r.x) : Math.ceil(r.x-1.0);
            r.y -= (r.y > 0.0) ? Math.floor(r.y) : Math.ceil(r.y-1.0);
 //           if(r.x > 0) {r.x -= Math.floor(r.x);}
 //           else         {r.x -= Math.ceil(r.x-1.0);}
 //           if(r.y > 0) {r.y -= Math.floor(r.y);}
 //           else         {r.y -= Math.ceil(r.y-1.0);}
        }
        public void inflate(double scale) {dimensions.TE(scale);}
        public double volume() {return dimensions.x * dimensions.y;}
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
    }  //end of BoundarySquarePeriodic
            
    private static final class AtomPair extends simulate.AtomPair {  
        AtomCoordinate c1;
        AtomCoordinate c2;
        final Boundary boundary;
        private final Vector dr = new Vector();
        private double drx, dry, dvx, dvy;
        public AtomPair() {boundary = new BoundaryNone();}
        public AtomPair(Boundary b) {boundary = b;}
        public AtomPair(Boundary b, Atom a1, Atom a2) {
            boundary = b;
            if(a1 != null && a2 != null) {
                c1 = (AtomCoordinate)a1.coordinate;  //cast from Space.AtomCoordinate
                c2 = (AtomCoordinate)a2.coordinate;
                reset();
            }
        }
        public void reset(Atom a1, Atom a2) {
            c1 = (AtomCoordinate)a1.coordinate;
            c2 = (AtomCoordinate)a2.coordinate;
            reset();
        }
        public void reset() {
            dr.x = c2.r.x - c1.r.x;
            dr.y = c2.r.y - c1.r.y;
//                dr.x -= (dr.x > 0.0) ? Math.floor(dr.x+0.5) : Math.ceil(dr.x-0.5);
//                dr.y -= (dr.y > 0.0) ? Math.floor(dr.y+0.5) : Math.ceil(dr.y-0.5);
            boundary.apply(dr);
            drx = dr.x; dry = dr.y;
            r2 = drx*drx + dry*dry;
            atom1 = c1.atom;
            atom2 = c2.atom;
            double rm1 = atom1.type.rm();
            double rm2 = atom2.type.rm();
            dvx = rm2*c2.p.x - rm1*c1.p.x;
            dvy = rm2*c2.p.y - rm1*c1.p.y;  
        }
                
//        public double r2() {
//            return drx*drx + dry*dry;
//        }
        public double v2() {
            double rm1 = c1.atom.type.rm();
            double rm2 = c2.atom.type.rm();
            return dvx*dvx + dvy*dvy;
        }
        public double vDotr() {
            double rm1 = c1.atom.type.rm();
            double rm2 = c2.atom.type.rm();
            return drx*dvx + dry*dvy;
        }
        public void push(double impulse) {  //changes momentum in the direction joining the atoms
            c1.p.x += impulse*drx;
            c1.p.y += impulse*dry;
            c2.p.x -= impulse*drx;
            c2.p.y -= impulse*dry;
        }
        public void setSeparation(double r2New) {
            double ratio = c2.atom.type.mass()*c1.atom.type.rm();  // (mass2/mass1)
            double delta = (Math.sqrt(r2New/r2()) - 1.0)/(1+ratio);
            c1.r.x -= ratio*delta*drx;
            c1.r.y -= ratio*delta*dry;
            c2.r.x += delta*drx;
            c2.r.y += delta*dry;
            //need call reset?
        }
   //     public final Atom atom1() {return c1.atom();}
   //     public final Atom atom2() {return c2.atom();}
    }

    public static final class Vector implements Space.Vector {  //declared final for efficient method calls
        public static final Random random = new Random();
        public static final Vector ORIGIN = new Vector(0.0,0.0);
        double x, y;
        public Vector () {x = 0.0; y = 0.0;}
        public Vector (double a1, double a2) {x = a1; y = a2;}
        public double component(int i) {return (i==0) ? x : y;}
        public void setComponent(int i, double d) {if(i==0) x=d; else y=d;}
        public void E(Space.Vector u) {E((Vector)u);}
        public void PE(Space.Vector u) {PE((Vector)u);}
        public void ME(Space.Vector u) {ME((Vector)u);}
        public void TE(Space.Vector u) {TE((Vector)u);}
        public void DE(Space.Vector u) {DE((Vector)u);}
        public double dot(Space.Vector u) {return dot((Vector)u);}
        public void E(Vector u) {x = u.x; y = u.y;}
        public void E(double a) {x = a; y = a;}
        public void Ea1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x = a1*u1.x; y = a1*u1.y;}
        public void PEa1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x += a1*u1.x; y += a1*u1.y;}
        public void PE(Vector u) {x += u.x; y += u.y;}
        public void ME(Vector u) {x -= u.x; y -= u.y;}
        public void TE(double a) {x *= a; y *= a;}
        public void DE(double a) {x /= a; y /= a;}
        public double squared() {return x*x + y*y;}
        public double dot(Vector u) {return x*u.x + y*u.y;}
        public void randomStep(double d) {x += (2.*random.nextDouble()-1.0)*d; y+= (2.*random.nextDouble()-1.0)*d;} //uniformly distributed random step in x and y, within +/- d
        public void setRandom(double d) {x = random.nextDouble()*d; y = random.nextDouble()*d;}
        public void setRandom(double dx, double dy) {x = random.nextDouble()*dx; y = random.nextDouble()*dy;}
        public void setRandom(Vector u) {setRandom(u.x,u.y);}
        public void setToOrigin() {x = ORIGIN.x; y = ORIGIN.y;}
        public void randomDirection() {x = Math.cos(2*Math.PI*random.nextDouble()); y = Math.sqrt(1.0 - x*x);}
    }
    
    static class Coordinate implements Space.Coordinate {
        public final Vector r = new Vector();  //Cartesian coordinates
        public final Vector p = new Vector();  //Momentum vector
        public Space.Vector position() {return r;}
        public Space.Vector momentum() {return p;}
        public double position(int i) {return r.component(i);}
        public double momentum(int i) {return p.component(i);}
    }    
    
    //much of AtomCoordinate is identical in every Space class
    //It is duplicated because it extends Coordinate, which is unique to each Space
    final static class AtomCoordinate extends Coordinate implements Space.AtomCoordinate {
        AtomCoordinate nextCoordinate, previousCoordinate;
        AtomCoordinate(Atom a) {atom = a;}  //constructor
        public final Atom atom;        
        private final Vector rLast = new Vector();
        private final Vector temp = new Vector();
        
        public void translateTo(Space.Vector u) {translateTo((Vector)u);}      
        public void translateBy(Space.Vector u) {translateBy((Vector)u);}
        public void translateToward(Space.Vector u, double d) {translateToward((Vector)u,d);}
        public void displaceTo(Space.Vector u) {displaceTo((Vector)u);}  //want to eliminate these casts
        public void displaceBy(Space.Vector u) {displaceBy((Vector)u);}
        public void accelerate(Space.Vector u) {accelerate((Vector)u);}
        public void accelerateToward(Space.Vector u, double d) {accelerateToward((Vector)u,d);}
        
        public void translateToward(Vector u, double d) {r.x += d*u.x; r.y += d*u.y;}
        public void displaceWithin(double d) {rLast.x = r.x; rLast.y = r.y; r.randomStep(d);}
        public void translateTo(Vector u) {r.x = u.x; r.y = u.y;}      
        public void translateBy(Vector u) {r.x += u.x; r.y += u.y;}  //no PBC is routinely applied
        public void displaceTo(Vector u) {rLast.x = r.x; rLast.y = r.y; r.x = u.x; r.y = u.y;}  
        public void displaceBy(Vector u) {rLast.x = r.x; rLast.y = r.y; r.x += u.x; r.y += u.y;}
        public void displaceToRandom(simulate.Phase p) {displaceTo(p.boundary().randomPosition());}
        public void translateToRandom(simulate.Phase p) {translateTo(p.boundary().randomPosition());}
        public void accelerateBy(Vector u) {p.x += u.x; p.y += u.y;}
        public void accelerateToward(Vector u, double d) {p.x += d*u.x; p.y += d*u.y;}
        public void replace() {r.x = rLast.x; r.y = rLast.y;}
        public void inflate(double s) {r.x *= s; r.y *= s;}
        public double kineticEnergy() {return 0.5*(p.x*p.x + p.y*p.y)*atom.type.rm();}
        public void randomizeMomentum(double temperature) {  //not very sophisticated; random only in direction, not magnitude
            double momentum = Math.sqrt(atom.mass()*temperature*(double)D/Constants.KE2T);  //need to divide by sqrt(m) to get velocity
            p.randomDirection();
            p.TE(momentum);
        }
        public void scaleMomentum(double scale) {p.x *= scale; p.y *= scale;}
        public Space.Vector velocity() {temp.E(p); temp.TE(atom.type.rm()); return temp;}  //returned vector is not thread-safe

/*        public void translateTo(Space.Vector u) {translateTo((Vector)u);}      
        public void translateBy(Space.Vector u) {r.PE((Vector)u);}
        public void translateToward(Space.Vector u, double d) {translateToward((Vector)u,d);}
        public void translateToward(Vector u, double d) {temp.Ea1Tv1(d,u); translateBy(temp);}
        public void displaceTo(Space.Vector u) {displaceTo((Vector)u);}  //want to eliminate these casts
        public void displaceBy(Space.Vector u) {displaceBy((Vector)u);}
        public void displaceWithin(double d) {temp.setToOrigin(); temp.randomStep(d); displaceBy(temp);}
//        public void displaceToRandom() {rLast.E(r); r.setRandom(dimensions.x, dimensions.y);} 
        public void accelerate(Space.Vector u) {p.PE((Vector)u);}
        public void translateTo(Vector u) {r.E(u);}      
        public void translateBy(Vector u) {r.PE(u);}  //no PBC is routinely applied
//        public void translateToRandom() {r.setRandom(dimensions.x, dimensions.y);}
        public void displaceTo(Vector u) {rLast.E(r); translateTo(u);}  
        public void displaceBy(Vector u) {rLast.E(r); translateBy(u);}
        public void displaceToRandom(simulate.Phase p) {displaceTo(p.boundary().randomPosition());}
        public void translateToRandom(simulate.Phase p) {translateTo(p.boundary().randomPosition());}
        public void accelerateBy(Vector u) {p.PE(u);}
        public void accelerateToward(Space.Vector u, double d) {accelerateToward((Vector)u,d);}
        public void accelerateToward(Vector u, double d) {temp.Ea1Tv1(d,u); accelerateBy(temp);}
        public void replace() {r.E(rLast);}
        public void inflate(double s) {r.TE(s);}
        public double kineticEnergy() {return 0.5*p.squared()*atom.type.rm();}
        public void randomizeMomentum(double temperature) {  //not very sophisticated; random only in direction, not magnitude
            double momentum = Math.sqrt(atom.mass()*temperature*(double)D/Constants.KE2T);  //need to divide by sqrt(m) to get velocity
            p.randomDirection();
            p.TE(momentum);
        }
        public void scaleMomentum(double scale) {p.TE(scale);}
        public Space.Vector position() {return r;}
        public Space.Vector momentum() {return p;}
        public double position(int i) {return r.component(i);}
        public double momentum(int i) {return p.component(i);}
        public Space.Vector velocity() {temp.E(p); temp.TE(atom.type.rm()); return temp;}  //returned vector is not thread-safe
      */  
        //following methods are same in all Space classes
        public Space.AtomCoordinate nextCoordinate() {return nextCoordinate;}
        public Space.AtomCoordinate previousCoordinate() {return previousCoordinate;}
        public final void setNextCoordinate(Space.AtomCoordinate c) {
           nextCoordinate = (AtomCoordinate)c;
           if(c != null) {((AtomCoordinate)c).previousCoordinate = this;}
        }
        public final void clearPreviousCoordinate() {previousCoordinate = null;}
        public final Atom previousAtom() {
            Space.AtomCoordinate c = atom.coordinate.previousCoordinate();
            return (c==null) ? null : c.atom();
        }
        public final Atom nextAtom() {
            Space.AtomCoordinate c = atom.coordinate.nextCoordinate();
            return (c==null) ? null : c.atom();
        }
        public final Atom atom() {return atom;}
    } 

    //These iterators are identical in every Space class; they are repeated in each
    //because they make direct use of the Coordinate type in the class; otherwise casting would be needed
    // Perhaps interitance would work, but haven't tried it
        
    //"Full" --> Each iteration of inner loop begins with same first atom
    private static final class PairIteratorFull implements simulate.AtomPair.Iterator.A {
        final AtomPair pair;
        AtomCoordinate outer, inner;
        private AtomCoordinate iFirst, iLast, oFirst, oLast;
        private boolean hasNext;
        public PairIteratorFull(Space.Boundary b) {  //null constructor
            pair = new AtomPair((Boundary)b);
            hasNext = false;
        }  
        public PairIteratorFull(Space.Boundary b, Atom iF, Atom iL, Atom oF, Atom oL) {  //constructor
            pair = new AtomPair((Boundary)b);
            reset(iF,iL,oF,oL);
        }
        public void reset(Atom iL, Atom oF, Atom oL) {reset(oF,iL,oF,oL);}  //take inner and outer first atoms as same
        public void reset(Atom iF, Atom iL, Atom oF, Atom oL) {
            if(iF == null || oF == null) {hasNext = false; return;}
            iFirst = (AtomCoordinate)iF.coordinate; 
            iLast =  (iL==null) ? null : (AtomCoordinate)iL.coordinate; 
            oFirst = (AtomCoordinate)oF.coordinate; 
            oLast =  (oL==null) ? null : (AtomCoordinate)oL.coordinate;
            inner = iFirst;
            outer = oFirst;
            hasNext = true;
        }
        public simulate.AtomPair next() {
            if(!hasNext) {return null;}
            pair.c1 = outer;   //c1 should always be outer
            pair.c2 = inner;
            pair.reset();
            if(inner == iLast) {                                     //end of inner loop
                if(outer == oLast) {hasNext = false;}                //all done
                else {outer = outer.nextCoordinate; inner = iFirst;} //advance outer, reset inner
            }
            else {inner = inner.nextCoordinate;}
            return pair;
        }
        public final void allDone() {hasNext = false;}   //for forcing iterator to indicate it has no more pairs
        public boolean hasNext() {return hasNext;}
        public void reset() {reset(iFirst.atom(), iLast.atom(), oFirst.atom(), oLast.atom());}
    }
        
    //"Half" --> Each iteration of inner loop begins with atom after outer loop atom
    private static final class PairIteratorHalf implements simulate.AtomPair.Iterator.A {
        final AtomPair pair;
        AtomCoordinate outer, inner;
        private AtomCoordinate iFirst, iLast, oFirst, oLast;
        private boolean hasNext;
        public PairIteratorHalf(Space.Boundary b) {
            pair = new AtomPair((Boundary)b);
            hasNext = false;
        }
        public PairIteratorHalf(Space.Boundary b, Atom iL, Atom oF, Atom oL) {  //constructor
            pair = new AtomPair((Boundary)b);
            reset(iL,oF,oL);
        }
        public void reset(Atom iF, Atom iL, Atom oF, Atom oL) {reset(iL,oF,oL);} //ignore first argument
        public void reset(Atom iL, Atom oF, Atom oL) {
            if(oF == null) {hasNext = false; return;}
            iLast =  (iL==null) ? null : (AtomCoordinate)iL.coordinate; 
            oFirst =  (AtomCoordinate)oF.coordinate; 
            oLast =  (iL==null) ? null : (AtomCoordinate)oL.coordinate;
            outer = oFirst;
            inner = outer.nextCoordinate;
            hasNext = (inner != null);
        }
        public simulate.AtomPair next() {
            pair.c1 = outer;   //c1 should always be outer
            pair.c2 = inner;
            pair.reset();
            if(inner == iLast) {                                     //end of inner loop
                if(outer == oLast) {hasNext = false;}                //all done
                else {outer = outer.nextCoordinate; inner = outer.nextCoordinate;} //advance outer, reset inner
            }
            else {inner = inner.nextCoordinate;}
            return pair;
        }
        public final void allDone() {hasNext = false;}   //for forcing iterator to indicate it has no more pairs
        public boolean hasNext() {return hasNext;}
        public void reset() {reset(iLast.atom(), oFirst.atom(), oLast.atom());}
    }    
}