package simulate;

public interface AtomPairIterator {
    
    public boolean hasNext();
    public AtomPair next();
    
    public interface SS extends AtomPairIterator {public void reset(Species s1, Species s2);}
    public interface S extends AtomPairIterator {public void reset(Species s);}
    public interface MM extends AtomPairIterator {public void reset(Molecule m1, Molecule m2);}
    public interface M extends AtomPairIterator {public void reset(Molecule m);}
    public interface A extends AtomPairIterator {
        public void reset(Atom a1, Atom a2, Atom a3);
        public void reset(Atom a1, Atom a2, Atom a3, Atom a4);
        public void allDone();
    }
}