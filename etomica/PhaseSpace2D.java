package simulate;

public class PhaseSpace2D extends PhaseSpace {
    
    public PhaseSpace2D() {
        dimensions = new PhaseSpace2D.Vector();
        dimensions.x = 1.0;
        dimensions.y = 1.0;
        computeVolume();
        periodic = false;
    }
 
    public PhaseSpace.AtomCoordinate makeAtomCoordinate(Atom a) {return new AtomCoordinate(a);}
    public PhaseSpace.MoleculeCoordinate makeMoleculeCoordinate(Molecule m) {return new MoleculeCoordinate(m);}
    public PhaseSpace.AtomPair makeAtomPair(Atom a1, Atom a2) {return new AtomPair(a1, a2);}
    
    public final AtomPairIterator.A makePairIteratorFull(Atom iF, Atom iL, Atom oF, Atom oL) {return new PairIteratorFull(iF,iL,oF,oL);}
    public final AtomPairIterator.A makePairIteratorHalf(Atom iL, Atom oF, Atom oL) {return new PairIteratorHalf(iL,oF,oL);}
    public final AtomPairIterator.A makePairIteratorFull() {return new PairIteratorFull();}
    public final AtomPairIterator.A makePairIteratorHalf() {return new PairIteratorHalf();}
 
    public static final class Vector implements PhaseSpace.Vector {
        double x, y;
        public Vector () {x = 0.0; y = 0.0;}
        public Vector (double a1, double a2) {x = a1; y = a2;}
        public void E(Vector u) {x = u.x; y = u.y;}
        public void E(double a) {x = a; y = a;}
        public void PE(Vector u) {x += u.x; y += u.y;}
        public void TE(double a) {x *= a; y *= a;}
        public void DE(double a) {x /= a; y /= a;}
        public double square() {return x*x + y*y;}
        public double dot(Vector u) {return x*u.x + y*u.y;}
    }
    
    abstract class Coordinate implements PhaseSpace.Coordinate {
        public final Vector r = new Vector();  //Cartesian coordinates
        public final Vector p = new Vector();  //Momentum vector
        public double mass;
    }    
    
    //much of AtomCoordinate and MoleculeCoordinate are identical in every PhaseSpace class
    //They are duplicated because they extend Coordinate, which is unique to each PhaseSpace
    final class AtomCoordinate extends Coordinate implements PhaseSpace.AtomCoordinate {
        AtomCoordinate nextCoordinate, previousCoordinate;
        AtomCoordinate(Atom a) {atom = a;}  //constructor
        public final Atom atom;
        
        protected final Vector rLast = new Vector();
        protected final Vector temp = new Vector();
        public void translateTo(Vector u) {r.E(u);}      //if using PBC, apply here
        public void translateBy(Vector u) {r.PE(u);}
        public void displaceTo(Vector u) {rLast.E(r); r.E(u);}
        public void displaceBy(Vector u) {rLast.E(r); r.PE(u);}
        public void replace() {r.E(rLast);}
        public void inflate(double s) {r.TE(s);}
        public void accelerate(Vector u) {p.PE(u);}
        public double kineticEnergy() {return 0.5*p.square()/mass;}
        public Vector position() {return r;}
        public Vector momentum() {return p;}
        public Vector velocity() {temp.E(p); temp.DE(mass); return temp;}  //returned vector is not thread-safe
        
        //following methods are same in all PhaseSpace classes
        public final void setNextCoordinate(PhaseSpace.Coordinate c) {
           nextCoordinate = (AtomCoordinate)c;
           if(c != null) {((AtomCoordinate)c).previousCoordinate = this;}
        }
        public final void clearPreviousCoordinate() {previousCoordinate = null;}
        public final Atom previousAtom() {
            AtomCoordinate c = atom.coordinate.previousCoordinate;
            return (c==null) ? null : c.atom;
        }
        public final Atom nextAtom() {
            AtomCoordinate c = atom.coordinate.nextCoordinate;
            return (c==null) ? null : c.atom;
        }
    }    
    final class MoleculeCoordinate extends Coordinate implements PhaseSpace.MoleculeCoordinate {
        MoleculeCoordinate nextCoordinate, previousCoordinate;
        MoleculeCoordinate(Molecule m) {molecule = m;}  //constructor
        public final Molecule molecule;
        public final void setNextCoordinate(PhaseSpace.Coordinate c) {
           nextCoordinate = (MoleculeCoordinate)c;
           if(c != null) {((MoleculeCoordinate)c).previousCoordinate = this;}
        }
        public final void clearPreviousCoordinate() {previousCoordinate = null;}
        public final Molecule previousMolecule() {
            MoleculeCoordinate c = molecule.coordinate.previousCoordinate;
            return (c==null) ? null : c.molecule;
        }
        public final Molecule nextMolecule() {
            MoleculeCoordinate c = molecule.coordinate.nextCoordinate;
            return (c==null) ? null : c.molecule;
        }
    }    
    
    private class AtomPair implements PhaseSpace.AtomPair {  //Inner AtomPair class
        Coordinate c1;
        Coordinate c2;
        public AtomPair() {}
        public AtomPair(Atom a1, Atom a2) {
            if(a1 != null && a2 != null) {
                c1 = (Coordinate)a1.coordinate;
                c2 = (Coordinate)a2.coordinate;
            }
        }
        
        public double r2() {
            double dx = c1.r.x - c2.r.x;   //change for PBC
            double dy = c1.r.y - c2.r.y;
            return dx*dx + dy*dy;
        }
        public final Atom atom1() {return c1.atom;}
        public final Atom atom2() {return c2.atom;}
    }
    
    //These iterators are identical in every PhaseSpace class; they are repeated in each
    //because they make direct use of the Coordinate type in the class; otherwise casting would be needed
    // Perhaps interitance would work, but haven't tried it
    
    //"Full" --> Each iteration of inner loop begins with same first atom
    private class PairIteratorFull implements AtomPairIterator.A {
        final AtomPair pair = new AtomPair();
        Coordinate outer, inner;
        private Coordinate iFirst, iLast, oLast;
        private boolean hasNext;
        public PairIteratorFull() {hasNext = false;}  //null constructor
        public PairIteratorFull(Atom iF, Atom iL, Atom oF, Atom oL) {reset(iF,iL,oF,oL);}  //constructor
        public void reset(Atom iL, Atom oF, Atom oL) {reset(oF,iL,oF,oL);}  //take inner and outer first atoms as same
        public void reset(Atom iF, Atom iL, Atom oF, Atom oL) {
            iFirst = (Coordinate)iF.coordinate; 
            iLast =  (Coordinate)iL.coordinate; 
            outer = (Coordinate)oF.coordinate; 
            oLast =  (Coordinate)oL.coordinate;
            inner = iFirst;
            hasNext = (inner != null && outer != null);
        }
        public AtomPair next() {
            if(!hasNext) {return null;}
            pair.c1 = outer;
            pair.c2 = inner;
            if(inner == iLast) {                                     //end of inner loop
                if(outer == oLast) {hasNext = false;}                //all done
                else {outer = outer.nextCoordinate; inner = iFirst;} //advance outer, reset inner
            }
            return pair;
        }
        public final void allDone() {hasNext = false;}   //for forcing iterator to indicate it has no more pairs
        public boolean hasNext() {return hasNext;}
    }
    
    //"Half" --> Each iteration of inner loop begins with atom after outer loop atom
    private class PairIteratorHalf implements AtomPairIterator.A {
        final AtomPair pair = new AtomPair();
        Coordinate outer, inner;
        private Coordinate iFirst, iLast, oLast;
        private boolean hasNext;
        public PairIteratorHalf() {hasNext = false;}
        public PairIteratorHalf(Atom iL, Atom oF, Atom oL) {reset(iL,oF,oL);}  //constructor
        public void reset(Atom iF, Atom iL, Atom oF, Atom oL) {reset(iL,oF,oL);} //ignore first argument
        public void reset(Atom iL, Atom oF, Atom oL) {
            if(!hasNext) {return null;}
            iLast =  (Coordinate)iL.coordinate; 
            outer = (Coordinate)oF.coordinate; 
            oLast =  (Coordinate)oL.coordinate;
            inner = outer.nextCoordinate;
            hasNext = (inner != null && outer != null);
        }
        public AtomPair next() {
            pair.c1 = outer;
            pair.c2 = inner;
            if(inner == iLast) {                                     //end of inner loop
                if(outer == oLast) {hasNext = false;}                //all done
                else {outer = outer.nextCoordinate; inner = outer.nextCoordinate;} //advance outer, reset inner
            }
            return pair;
        }
        public final void allDone() {hasNext = false;}   //for forcing iterator to indicate it has no more pairs
        public boolean hasNext() {return hasNext;}
    }    

 /**
  * Size of Phase (width, height) in Angstroms
  * Default value is 1.0 for each dimension.
  */
    private final Vector dimensions = new Vector(1.0,1.0);


}