package simulate;

public interface AtomPairIterator {
    public boolean hasNext();
//    public boolean hasPrevious();
//    public boolean hasNextInSpecies(Species s);
//    public boolean hasPreviousInSpecies(Species s);
//    public boolean hasNextInMolecule(Molecule m);
//    public boolean hasPreviousInMolecule(Molecule m);
    public SpaceAtomPair next();
//    public SpaceAtomPair previous();
//    public SpaceAtomPair nextInMolecule();
//    public void reset();
//    public void reset(SpaceAtom a1, SpaceAtom a2);  //sets to first pair in species
//    public void reset(Molecule m);
//    public void reset(Species s);
//    public void reset(Species s1, Species s2);
}