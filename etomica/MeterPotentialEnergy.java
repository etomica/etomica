package simulate;

public class MeterPotentialEnergy extends simulate.Meter
{
  AtomPair.Iterator iterator;
  AtomPair.Iterator.MP iteratorMP;
  AtomPair.Iterator.All iteratorAll;
  
  AtomPair.Iterator apiUp, apiDown;
  
  public MeterPotentialEnergy() {
    super();
    setLabel("Potential Energy");
  }
  
  public void setPhase(Phase p) {
    super.setPhase(p);
    iteratorMP = new AtomPair.Iterator.MP(p);
    iteratorAll = p.iterator.makeAtomPairIteratorAll();
    
    apiUp = p.iterator.makeAtomPairIteratorUp();
    apiDown = p.iterator.makeAtomPairIteratorDown();
  }
    

 /**
  * Computes total potential energy for all atom pairs in phase
  * Returns infinity (MAX_VALUE) as soon as overlap is detected
  */
    public final double currentValue() {
                        //could make special case for single species, monatomic
        double pe = 0.0;
        iteratorAll.reset();
        while(iteratorAll.hasNext()) {
            AtomPair pair = iteratorAll.next();
            double energy = phase.parentSimulation.getPotential(pair).energy(pair);
            if(energy == Double.MAX_VALUE) return Double.MAX_VALUE;
            pe += energy;
        }
        return pe;
    }


  /**
   * Computes intermolecular contribution to potential energy for atom.  Returns zero if atom is not in phase.
   * Does not include intramolecular contribution
   * Returns infinity as soon as overlap is encountered

   //  May want to change this so it handles atoms not in phase

   */
    public final double currentValue(Atom a) {
        if(phase != a.parentPhase()) {return 0.0;}  //also handles condition that phase contains no atoms
        double pe = 0.0;
        apiUp.reset(a);
        while(apiUp.hasNext()) {
            AtomPair pair = apiUp.next();
            double energy = phase.parentSimulation.getPotential(pair).energy(pair);
            if(energy == Double.MAX_VALUE) return Double.MAX_VALUE;
            pe += energy;
        }
        apiDown.reset(a);
        while(apiDown.hasNext()) {
            AtomPair pair = apiDown.next();
            double energy = phase.parentSimulation.getPotential(pair).energy(pair);
            if(energy == Double.MAX_VALUE) return Double.MAX_VALUE;
            pe += energy;
        }
        return pe;
    }


 /**
  * Computes and returns the intermolecular contribution of a molecule to energy of phase
  * Does not include intramolecular energy of molecule
  * @return  this molecule's inter-molecular potential energy, divided by kB, in Kelvin
  */
    
    public final double currentValue(Molecule m) {
        iteratorMP.reset(m);
        double pe = 0.0;
        while(iteratorMP.hasNext()) {            //iterator not returning a pair here
            AtomPair pair = iteratorMP.next();
//            pe += pair.potential.energy(pair);
            pe += phase.parentSimulation.getPotential(pair).energy(pair);
        }
        return pe;
    }
         //looks like these methods are exactly the same************
 /**
  * Computes and returns the total intermolecular energy of a molecule with all molecules in phase
  * Does not include intramolecular energy of molecule
  */
    
    public final double insertionValue(Molecule m) {

        iteratorMP.reset(m);
        double pe = 0.0;
        while(iteratorMP.hasNext()) {
            AtomPair pair = iteratorMP.next();
            double energy = phase.parentSimulation.getPotential(pair).energy(pair);
            if(energy == Double.MAX_VALUE) return Double.MAX_VALUE;
            pe += energy;
        }    
      return pe;
/*      AtomPair p = new AtomPair(phase);
      int count = 0;
      for(Atom a=phase.lastAtom(); a!=null; a=a.previousAtom()) {
            count++;
            p.reset(m.firstAtom,a);
            double energy = phase.parentSimulation.getPotential(p).energy(p);
            if(energy == Double.MAX_VALUE) {
                System.out.println(count+" "+phase.atomCount);
                return Double.MAX_VALUE;
            }
            pe += energy;
      }
      return pe; */  //debugging --- delete
    }

    
}
