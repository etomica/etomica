package simulate;

public abstract class Space {
    
    public abstract int D();
    
    public abstract Vector makeVector();      //Space.Vector
    public abstract Coordinate makeCoordinate(Occupant o);
    public abstract CoordinatePair makeCoordinatePair(Coordinate c1, Coordinate c2, Boundary b);
    public abstract CoordinatePair makeCoordinatePair(Boundary b);
    public abstract Boundary makeBoundary(int iBoundary);

    interface Occupant {
        public Coordinate coordinate();
        public double mass();
        public double rm();
    }
    
//  Vector contains what is needed to describe a point in the space
    interface Vector {    //probably want this to be an abstract class
        public double component(int i);
        public void setComponent(int i, double d);
        public void E(Vector u);
        public void E(double a);
        public void PE(Vector u);
        public void ME(Vector u);
        public void TE(Vector u);
        public void DE(Vector u);
        public void TE(double a);
        public void DE(double a);
        public void Ea1Tv1(double a, Vector u);
        public double squared();
        public double dot(Vector u);
        public void setRandom(double d);
        public void PEa1Tv1(double a, Vector u);
        public void setRandomSphere();  //random point in unit sphere
        public void setRandomCube(); //random point in a unit cube
    }

//  Coordinate collects all vectors needed to describe point in phase space -- position and (maybe) momentum
    public static abstract class Coordinate {
        protected final Space.Occupant parent;        
        Coordinate(Occupant p) {parent = p;}          //constructor
        public final Space.Occupant parent() {return parent;}
        public abstract Vector makeVector();
        public abstract Vector position();
        public abstract Vector momentum();
        public abstract double position(int i);
        public abstract double momentum(int i);
        public abstract double kineticEnergy(double mass);
        public void scaleMomentum(double scale) {momentum().TE(scale);}
/*        public void setNextCoordinate(Coordinate c);
        public void clearPreviousCoordinate();
        
        public void translateTo(Vector r);
        public void translateToward(Vector e, double amount);
        public void translateBy(Vector dr);
        public void displaceTo(Vector r);
//        public void displaceBy(Vector dr);
        public void displaceWithin(double d);
        public void displaceToRandom(Phase p);
        public void translateToRandom(Phase p);
        public void randomizeMomentum(double temperature);
        public void replace();
        public void inflate(double s);
        public double kineticEnergy();
        public Vector position();
        public Vector momentum();

atomCoordinate methods
        public void translateTo(Vector r);
        public void translateToward(Vector e, double amount);
        public void translateBy(Vector dr);
        public void displaceTo(Vector r);
        public void displaceBy(Vector dr);
        public void displaceWithin(double d);
        public void displaceToRandom(Phase p);
        public void translateToRandom(Phase p);
        public void randomizeMomentum(double temperature);
        public void replace();
        public void inflate(double s);
        public double kineticEnergy();
        public Vector position();
        public Vector momentum();
        public double position(int i);
        public double momentum(int i);
 
        public Atom nextAtom();
        public Atom previousAtom();
        public Atom atom();
        public AtomCoordinate nextCoordinate();
        public AtomCoordinate previousCoordinate();
        public void scaleMomentum(double scale);
        public void accelerate(Vector dv);
        public void accelerateToward(Vector e, double amount);
        public Vector velocity();
        public void setNextCoordinate(AtomCoordinate c);
        public void clearPreviousCoordinate();
        public void setNextNeighbor(Space.AtomCoordinate c);
        public void clearPreviousNeighbor();
        public Space.AtomCoordinate nextNeighbor();
        public Space.AtomCoordinate previousNeighbor();
        public void assignCell();
        */
    }
    
    public static abstract class CoordinatePair {
        public Coordinate coordinate1, coordinate2;
        public double r2;
        public Potential potential;
        public CoordinatePair() {}  //null constructor
//        public CoordinatePair(Boundary b) {boundary = b;}
//        public CoordinatePair(Space.Occupant o1, Space.Occupant o2, Boundary b)
        public abstract void reset();
        public abstract void reset(Space.Coordinate c1, Space.Coordinate c2);
        public abstract double v2();
        public abstract double vDotr();
        public abstract void push(double impulse);  //impart equal and opposite impulse to momenta
        public abstract void setSeparation(double r2New);  //set square-distance between pair to r2New, by moving them along line joining them, keeping center of mass unchanged
        public final Coordinate coordinate1() {return coordinate1;}
        public final Coordinate coordinate2() {return coordinate2;}
        public final double r2() {return r2;}       
    }

    interface NeighborIterator {
        public void setNeighborRadius(double radius);
        public double getNeighborRadius();
        public Iterator makeIterator(Phase p);
        public void clear();
    }
        
    interface Boundary {
        public static final int NONE = 0;
        public static final int PERIODIC = 1;
        public static final int DEFAULT = PERIODIC; //default PBC is periodic
        public void centralImage(Vector r);
        public double volume();
        public Vector dimensions();
        public Vector randomPosition();
        public double[][] getOverflowShifts(Vector r, double distance);
        public void inflate(double s);
    /** Set of vectors describing the displacements needed to translate the central image
        *  to all of the periodic images.  Returns a two dimensional array of doubles.  The
        *  first index specifies each perioidic image, while the second index indicates the
        *  x and y components of the translation vector.
        *  Likely to override in subclasses.
        *
        *  @param nShells the number of shells of images to be computed
        */
        public abstract double[][] imageOrigins(int nShells);
    }
    
    // No-op version of boundary for manipulation at design time
    // Used by Phase.Virtual
//    public static class VirtualBoundary implements Boundary {
//        public void centralImage(Vector r) {}
//        public double volume() {return 0.0;}
//        public Vector dimensions() {return null;}
//        public Vector randomPosition() {return null;}
//       public double[][] getOverflowShifts(Vector r, double distance) {return null;}
//        public void inflate(double s) {}
//        public  double[][] imageOrigins(int nShells) {return null;}
//    }

    public Potential makePotential() {return new PotentialIdealGas();}  //default  
    
    
}    