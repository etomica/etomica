package simulate;

public class MeterPotentialEnergy extends simulate.Meter
{
  SpaceAtomPairIterator iterator;
  
  public MeterPotentialEnergy() {
    super();
    setLabel("Potential Energy");
  }

 /**
  * Computes total potential energy for all atom pairs in phase
  */
    public final double currentValue() {  
      double pe = 0.0;
      for(Species s1=phase.firstSpecies(); s1!=null; s1=s1.getNextSpecies()) {
        Potential1 p1 = phase.potential1[s1.speciesIndex];
        allIntra.reset(s1);
        while(allIntra.hasNext()) {
            SpaceAtomPair pair = allIntra.next();
            pe += p1.getPotential(pair.atom1(),pair.atom2()).energy(pair);
        }
        for(Species s2=s1; s2!=null; s2=s2.getNextSpecies()) {
            Potential2 p2 = phase.potential2[s1.speciesIndex][s2.speciesIndex];
            allInter.reset(s1,s2);
            while(allInter.hasNext()) {
                SpaceAtomPair pair = allInter.next();
                pe += p2.getPotential(pair.atom1(),pair.atom2()).energy(pair);
            }
        }
      }
      return pe;
    }
  
  /**
   * Computes total potential energy for atom, which must be in the list of atoms in phase
   */
    public final double currentValue(AtomC a) {
      double pe = 0.0;
      if(a.parentMolecule.nAtoms > 1) {pe += upListIntra(a) + downListIntra(a);}
      pe += upListInter(a) + downListInter(a);
      return pe;
    }
    
  /**
   * Computes total potential energy for atom (not present in phase), with all atoms in phase
   * Does not check to see if atom is actually in phase, which would contribute a (probably catastrophic) self-ineraction
   * Does not include intramolecular potential of atom
   */
    public final double insertionValue(AtomC a) {
      double pe = 0.0;
      Potential2[] p2 = phase.potential2[a.getSpeciesIndex()];
      for(AtomC b=(AtomC)phase.firstAtom(); b!=null; b=b.getNextAtomC()) {       //intermolecular
          if(b == a) System.out.println("atom error");
          pe += p2[b.getSpeciesIndex()].getPotential(a,b).energy(a,b);
      }
      return pe;
    }
    
 /**
  * Computes and returns the total (intra- and inter-molecular) contribution of a molecule
  * Includes intramolecular energy of molecule
  * @return  this molecule's intra- and inter-molecular potential energy, divided by kB, in Kelvin
  */
    public final double currentValue(Molecule m) {
      double pe = 0.0;
      Atom nextMoleculeAtom = m.terminationAtom();
      for(AtomC a=(AtomC)m.firstAtom(); a!=nextMoleculeAtom; a=a.getNextAtomC()) {  //loop over atoms in molecule
          pe += upListIntra(a) + upListInter(a) + downListInter(a);
      }
      return pe;
    }
    
 /**
  * Computes and returns the total intermolecular energy of a molecule (not present in phase) with all molecules in phase
  * Does not check to see if molecule is actually in phase, which would contribute a (probably catastrophic) self-ineraction
  * Does not include intramolecular energy of molecule
  */
    public final double insertionValue(Molecule m) {
      double pe = 0.0;
      Atom nextMoleculeAtom = m.terminationAtom();
      for(AtomC a=(AtomC)m.firstAtom(); a!=nextMoleculeAtom; a=a.getNextAtomC()) {  //loop over atoms in molecule
          pe += insertionValue(a);
      }
      return pe;
    }

    
    private double upListInter(AtomC a) {
      double pe = 0.0;
      AtomC nextMoleculeAtom = (AtomC)a.nextMoleculeFirstAtom();  //first atom on next molecule
      Potential2[] p2 = phase.potential2[a.getSpeciesIndex()];
      for(AtomC b=nextMoleculeAtom; b!=null; b=b.getNextAtomC()) {       //intermolecular
          pe += p2[b.getSpeciesIndex()].getPotential(a,b).energy(a,b);
      }
      return pe;
    }

    private double downListInter(AtomC a) {
      double pe = 0.0;
      AtomC previousMoleculeAtom = (AtomC)a.previousMoleculeLastAtom();  //first atom on next molecule
      Potential2[] p2 = phase.potential2[a.getSpeciesIndex()];
      for(AtomC b=previousMoleculeAtom; b!=null; b=b.getPreviousAtomC()) {       //intermolecular
          pe += p2[b.getSpeciesIndex()].getPotential(a,b).energy(a,b);
      }
      return pe;
    }

    private double upListIntra(AtomC a) {
      double pe = 0.0;
      AtomC nextMoleculeAtom = (AtomC)a.nextMoleculeFirstAtom();  //first atom on next molecule
      Potential1 p1 = phase.potential1[a.getSpeciesIndex()];            
      for(AtomC b=a.getNextAtomC(); b!=nextMoleculeAtom; b=b.getNextAtomC()) { //intramolecular
          pe += p1.getPotential(a,b).energy(a,b);
      }
      return pe;
    }

    private double downListIntra(AtomC a) {
      double pe = 0.0;
      AtomC previousMoleculeAtom = (AtomC)a.previousMoleculeLastAtom();  //first atom on next molecule
      Potential1 p1 = phase.potential1[a.getSpeciesIndex()];            
      for(AtomC b=a.getPreviousAtomC(); b!=previousMoleculeAtom; b=b.getPreviousAtomC()) { //intramolecular
          pe += p1.getPotential(a,b).energy(a,b);
      }
      return pe;
    }

}
