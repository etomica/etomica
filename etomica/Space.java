package simulate;

public abstract class Space {
    
    public abstract int D();
    
    public abstract Space.AtomCoordinate makeAtomCoordinate(Atom a);      //Space prefix is redundant
    public abstract Space.MoleculeCoordinate makeMoleculeCoordinate(Molecule m);
    public abstract Space.Vector makeVector();
    public abstract Phase makePhase(int boundary);

//  Vector contains what is needed to describe a point in the space
    interface Vector {
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
        public double squared();
        public double dot(Vector u);
    }

//  Coordinate collects all vectors needed to describe point in phase space -- position and (maybe) momentum
    interface Coordinate {
        public void setNextCoordinate(Coordinate c);
        public void clearPreviousCoordinate();
        
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
    }
    
    interface AtomCoordinate extends Coordinate {      //cannot be a class here because must inherit from Coordinate as it is defined in the PhaseSpace subclass
        public Atom nextAtom();
        public Atom previousAtom();
        public Atom atom();
        public AtomCoordinate nextCoordinate();
        public AtomCoordinate previousCoordinate();
        public void scaleMomentum(double scale);
        public void accelerate(Vector dv);
        public void accelerateToward(Vector e, double amount);
        public Vector velocity();
    }
    interface MoleculeCoordinate extends Coordinate {
        public Molecule nextMolecule();
        public Molecule previousMolecule();
        public Molecule molecule();
        public MoleculeCoordinate nextCoordinate();
        public MoleculeCoordinate previousCoordinate();
    }
    interface Boundary {
        public void centralImage(Vector r);
        public double volume();
        public Vector dimensions();
        public Vector randomPosition();
        public double[][] getOverflowShifts(Vector r, double distance);
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
}    