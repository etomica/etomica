package simulate;

public class PhaseSpace2D extends PhaseSpace {
    
    public PhaseSpace2D() {}
 
    public PhaseSpace.Coordinate makeCoordinate(Atom atom) {
        return new Coordinate(atom);
    }
    
    public final AtomPairIterator.A makePairIteratorFull(Atom iF, Atom iL, Atom oF, Atom oL) {return new PairIteratorFull(iF,iL,oF,oL);}
    public final AtomPairIterator.A makePairIteratorHalf(Atom iL, Atom oF, Atom oL) {return new PairIteratorHalf(iL,oF,oL);}
    public final AtomPairIterator.A makePairIteratorFull() {return new PairIteratorFull();}
    public final AtomPairIterator.A makePairIteratorHalf() {return new PairIteratorHalf();}
 
    class Coordinate implements PhaseSpace.Coordinate {
        Coordinate(Atom a) {atom = a;}  //constructor
        public final Atom atom;
        public final Vector r = new Vector();  //Cartesian coordinates
        public final Vector p = new Vector();  //Momentum vector
        public final SpaceVector makeVector() {return new Vector();}
        Coordinate nextCoordinate;
        Coordinate previousCoordinate;
        public final void setNextCoordinate(PhaseSpaceCoordinate c) {
           nextCoordinate = (Coordinate)c;
           if(c != null) {((Coordinate)c).previousCoordinate = this;}
        }
        public final void clearPreviousCoordinate() {previousCoordinate = null;}
    }    
    
    public class Vector implements SpaceVector {
        double x = 0.0;
        double y = 0.0; 
    }
    
    private class IAtomPair implements AtomPair {  //Inner AtomPair class
        Coordinate c1;
        Coordinate c2;
        
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
        final IAtomPair pair = new IAtomPair();
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
        final IAtomPair pair = new IAtomPair();
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
}