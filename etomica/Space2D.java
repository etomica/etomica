package simulate;
import java.awt.Graphics;
import java.util.Random;

public class Space2D extends Space {
    
    public static final int D = 2;
    public final int D() {return D;}

    public Space.AtomCoordinate makeAtomCoordinate(Atom a) {return new AtomCoordinate(a);}
    public Space.MoleculeCoordinate makeMoleculeCoordinate(Molecule m) {return new MoleculeCoordinate(m);}
    public Space.Vector makeVector() {return new Vector();}
    public simulate.Phase makePhase(int b) {
        switch(b) {
            case(Boundary.NONE):            return new Phase(new BoundaryNone());
            case(Boundary.PERIODIC_SQUARE): return new Phase(new BoundaryPeriodicSquare());
            default:                        return null;
        }
    }
    
    public static class Phase extends simulate.Phase {
        private Boundary boundary;
        public Phase(Boundary b) {
            boundary = b;
        }
        public final Space.Boundary boundary() {return boundary;}
        public final double volume() {return boundary.volume();}  //infinite volume unless using PBC
        public void inflate(double scale) {boundary.inflate(scale);}
        public final Space.Vector dimensions() {return boundary.dimensions;}
        public final simulate.AtomPair.Iterator.A makePairIteratorFull(Atom iF, Atom iL, Atom oF, Atom oL) {return new PairIteratorFull(boundary,iF,iL,oF,oL);}
        public final simulate.AtomPair.Iterator.A makePairIteratorHalf(Atom iL, Atom oF, Atom oL) {return new PairIteratorHalf(boundary,iL,oF,oL);}
        public final simulate.AtomPair.Iterator.A makePairIteratorFull() {return new PairIteratorFull(boundary);}
        public final simulate.AtomPair.Iterator.A makePairIteratorHalf() {return new PairIteratorHalf(boundary);}

        public simulate.AtomPair makeAtomPair(Atom a1, Atom a2) {return new AtomPair(boundary, a1, a2);}
        public void paint(Graphics g, int[] origin, double scale) {}
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
            if(dr.x > 0) {dr.x -= (dr.x > +0.5) ? dimensions.x : 0.0;}
            else         {dr.x += (dr.x < -0.5) ? dimensions.x : 0.0;}
            if(dr.y > 0) {dr.y -= (dr.y > +0.5) ? dimensions.y : 0.0;}
            else         {dr.y += (dr.y < -0.5) ? dimensions.y : 0.0;}
            dr.x *= dimensions.x;
            dr.y *= dimensions.y;
        }
        public void centralImage(Vector r) {
            //central-image code
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
            
    private static final class AtomPair implements simulate.AtomPair {  
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
        public void reset() {
            dr.x = c2.r.x - c1.r.x;
            dr.y = c2.r.y - c1.r.y;
            boundary.apply(dr);
            drx = dr.x; dry = dr.y;
        }
                
        public double r2() {
            return drx*drx + dry*dry;
        }
        public double v2() {
            double rm1 = c1.atom.type.rm();
            double rm2 = c2.atom.type.rm();
            double dvx = rm2*c2.p.x - rm1*c1.p.x;
            double dvy = rm2*c2.p.y - rm1*c1.p.y;
            return dvx*dvx + dvy*dvy;
        }
        public double vDotr() {
            double rm1 = c1.atom.type.rm();
            double rm2 = c2.atom.type.rm();
            double dvx = rm2*c2.p.x - rm1*c1.p.x;
            double dvy = rm2*c2.p.y - rm1*c1.p.y;
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
        public final Atom atom1() {return c1.atom();}
        public final Atom atom2() {return c2.atom();}
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
        public void Ea1Tv1(double a1, Vector u1) {x = a1*u1.x; y = a1*u1.y;}
        public void PEa1Tv1(double a1, Vector u1) {x += a1*u1.x; y += a1*u1.y;}
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
    
    public static final class VectorInt {  //declared final for efficient method calls
        int x, y;
        public VectorInt () {x = 0; y = 0;}
        public VectorInt (int a1, int a2) {x = a1; y = a2;}
        public void E(VectorInt u) {x = u.x; y = u.y;}
        public void E(int a) {x = a; y = a;}
        public void PE(VectorInt u) {x += u.x; y += u.y;}
        public void TE(int a) {x *= a; y *= a;}
        public void DE(int a) {x /= a; y /= a;}
        public int squared() {return (x*x + y*y);}
        public int dot(VectorInt u) {return x*u.x + y*u.y;}
    }
    
    static abstract class Coordinate implements Space.Coordinate {
        public final Vector r = new Vector();  //Cartesian coordinates
        public final Vector p = new Vector();  //Momentum vector
    }    
    
    //much of AtomCoordinate and MoleculeCoordinate are identical in every Space class
    //They are duplicated because they extend Coordinate, which is unique to each Space
    final static class AtomCoordinate extends Coordinate implements Space.AtomCoordinate {
        AtomCoordinate nextCoordinate, previousCoordinate;
        AtomCoordinate(Atom a) {atom = a;}  //constructor
        public final Atom atom;        
        private final Vector rLast = new Vector();
        private final Vector temp = new Vector();
        
        public void translateTo(Space.Vector u) {translateTo((Vector)u);}      
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
        
        //following methods are same in all Space classes
        public Space.AtomCoordinate nextCoordinate() {return nextCoordinate;}
        public Space.AtomCoordinate previousCoordinate() {return previousCoordinate;}
        public final void setNextCoordinate(Space.Coordinate c) {
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
    
    final static class MoleculeCoordinate extends Coordinate implements Space.MoleculeCoordinate {
        MoleculeCoordinate nextCoordinate, previousCoordinate;
        MoleculeCoordinate(Molecule m) {molecule = m;}  //constructor
        public final Molecule molecule;

        protected final Vector rLast = new Vector();
        protected final Vector temp = new Vector();
        public void updateR() {  //recomputes COM position from atom positions
            AtomCoordinate c = (AtomCoordinate)molecule.firstAtom.coordinate;
            if(molecule.nAtoms==1) {r.E(c.r);}  //one atom in molecule
            else {  //multiatomic
                r.Ea1Tv1(c.atom.mass(),c.r);
                do {c=c.nextCoordinate; r.PEa1Tv1(c.atom.mass(),c.r);} while (c.atom!=molecule.lastAtom);
                r.DE(molecule.mass());
            }
        }
        public void updateP() {  //recomputes total momentum from atom momenta
            AtomCoordinate c = (AtomCoordinate)molecule.firstAtom.coordinate;
            p.E(c.p);
            if(molecule.nAtoms==1) {return;}  //one atom in molecule
            do {c=c.nextCoordinate; p.PE(c.p);} while (c.atom!=molecule.lastAtom);
        }
        public void translateToward(Space.Vector u, double d) {temp.Ea1Tv1(d,(Vector)u); translateBy(temp);}
        public void translateBy(Space.Vector u) {translateBy((Vector)u);}
        public void translateBy(Vector u) {
            AtomCoordinate c = (AtomCoordinate)molecule.firstAtom.coordinate;
            c.translateBy(u);
            if(molecule.nAtoms == 1) {return;}
            do {c=c.nextCoordinate; c.translateBy(u);} while (c.atom!=molecule.lastAtom);
        }
//        public void displaceToRandom(Space.Vector dim) {temp.setRandom((Vector)dim); displaceTo(temp);} 
        public void translateTo(Space.Vector u) {
            updateR();  //update COM vector
            temp.E((Vector)u);  //temp = destination vector
            temp.ME(r);   //temp = destination - original = dr
            translateBy(temp);
        }
        public void displaceBy(Space.Vector u) {displaceBy((Vector)u);}
        public void displaceBy(Vector u) {
            AtomCoordinate c = (AtomCoordinate)molecule.firstAtom.coordinate;
            c.displaceBy(u);
            if(molecule.nAtoms == 1) {return;}
            do {c=c.nextCoordinate; c.displaceBy(u);} while (c.atom!=molecule.lastAtom);
        }
        public void displaceTo(Space.Vector u) {displaceTo((Vector)u);}
        public void displaceTo(Vector u) {
            updateR();  //update COM vector
            temp.E(u);  //temp = destination vector
            temp.ME(r);   //temp = destination - original = dr
            displaceBy(temp);
        }
        public void displaceWithin(double d) {
            temp.setRandom(d);
            displaceBy(temp);
        }
            
        public void displaceToRandom(simulate.Phase p) {displaceTo(p.boundary().randomPosition());}
        public void translateToRandom(simulate.Phase p) {translateTo(p.boundary().randomPosition());}

        public void replace() {
            AtomCoordinate c = (AtomCoordinate)molecule.firstAtom.coordinate;
            c.replace();
            if(molecule.nAtoms == 1) {return;}
            do {c=c.nextCoordinate; c.replace();} while (c.atom!=molecule.lastAtom);
        }
        public void inflate(double s) {
            updateR();
            temp.Ea1Tv1(s-1.0,r);
            displaceBy(temp);   //displaceBy doesn't use temp
        }
        public Space.Vector position() {updateR(); return r;}
        public Space.Vector momentum() {updateP(); return p;}
        public double position(int i) {updateR(); return r.component(i);}
        public double momentum(int i) {updateR(); return r.component(i);}
        public double kineticEnergy() {return 0.5*p.squared()*molecule.rm();}
        public void randomizeMomentum(double temperature) {
            AtomCoordinate c = (AtomCoordinate)molecule.firstAtom.coordinate;
            c.randomizeMomentum(temperature);
            if(molecule.nAtoms == 1) {return;}
            do {c=c.nextCoordinate; c.randomizeMomentum(temperature);} while (c.atom!=molecule.lastAtom);
        }
        public Space.MoleculeCoordinate nextCoordinate() {return nextCoordinate;}
        public Space.MoleculeCoordinate previousCoordinate() {return previousCoordinate;}
        public final void setNextCoordinate(Space.Coordinate c) {
           nextCoordinate = (MoleculeCoordinate)c;
           if(c != null) {((MoleculeCoordinate)c).previousCoordinate = this;}
        }
        public final void clearPreviousCoordinate() {previousCoordinate = null;}
        public final Molecule previousMolecule() {
            Space.MoleculeCoordinate c = molecule.coordinate.previousCoordinate();
            return (c==null) ? null : c.molecule();
        }
        public final Molecule nextMolecule() {
            Space.MoleculeCoordinate c = molecule.coordinate.nextCoordinate();
            return (c==null) ? null : c.molecule();
        }
        public final Molecule molecule() {return molecule;}
    }    
    
    //These iterators are identical in every Space class; they are repeated in each
    //because they make direct use of the Coordinate type in the class; otherwise casting would be needed
    // Perhaps interitance would work, but haven't tried it
        
    //"Full" --> Each iteration of inner loop begins with same first atom
    private static class PairIteratorFull implements simulate.AtomPair.Iterator.A {
        final AtomPair pair;
        AtomCoordinate outer, inner;
        private AtomCoordinate iFirst, iLast, oLast;
        private boolean hasNext;
        public PairIteratorFull(Boundary b) {  //null constructor
            pair = new AtomPair(b);
            hasNext = false;
        }  
        public PairIteratorFull(Boundary b, Atom iF, Atom iL, Atom oF, Atom oL) {  //constructor
            pair = new AtomPair(b);
            reset(iF,iL,oF,oL);
        }
        public void reset(Atom iL, Atom oF, Atom oL) {reset(oF,iL,oF,oL);}  //take inner and outer first atoms as same
        public void reset(Atom iF, Atom iL, Atom oF, Atom oL) {
            if(iF == null || oF == null) {hasNext = false; return;}
            iFirst = (AtomCoordinate)iF.coordinate; 
            iLast =  (iL==null) ? null : (AtomCoordinate)iL.coordinate; 
            outer = (AtomCoordinate)oF.coordinate; 
            oLast =  (oL==null) ? null : (AtomCoordinate)oL.coordinate;
            inner = iFirst;
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
    }
        
    //"Half" --> Each iteration of inner loop begins with atom after outer loop atom
    private static class PairIteratorHalf implements simulate.AtomPair.Iterator.A {
        final AtomPair pair;
        AtomCoordinate outer, inner;
        private AtomCoordinate iFirst, iLast, oLast;
        private boolean hasNext;
        public PairIteratorHalf(Boundary b) {
            pair = new AtomPair(b);
            hasNext = false;
        }
        public PairIteratorHalf(Boundary b, Atom iL, Atom oF, Atom oL) {  //constructor
            pair = new AtomPair(b);
            reset(iL,oF,oL);
        }
        public void reset(Atom iF, Atom iL, Atom oF, Atom oL) {reset(iL,oF,oL);} //ignore first argument
        public void reset(Atom iL, Atom oF, Atom oL) {
            if(oF == null) {hasNext = false; return;}
            iLast =  (iL==null) ? null : (AtomCoordinate)iL.coordinate; 
            outer =  (AtomCoordinate)oF.coordinate; 
            oLast =  (iL==null) ? null : (AtomCoordinate)oL.coordinate;
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
    }    
}