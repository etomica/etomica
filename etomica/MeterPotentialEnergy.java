package simulate;

public class MeterPotentialEnergy extends simulate.Meter
{
  AtomPair.Iterator iterator;
  AtomPair.Iterator.MP iteratorMP;
  AtomPair.Iterator iteratorAll;
  
  AtomPair.Iterator.A apiUp, apiDown;
  
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
  */
    public final double currentValue() {
                        //could make special case for single species, monatomic
        double pe = 0.0;
        iteratorAll.reset();
        while(iteratorAll.hasNext()) {
            AtomPair pair = iteratorAll.next();
            pe += phase.parentSimulation.getPotential(pair).energy(pair);
        }
        return pe;
    }

  /**
   * Computes intermolecular contribution to potential energy for atom.  Returns zero if atom is not in phase.
   * Does not include intramolecular contribution

   //  May want to change this so it handles atoms not in phase

   */
    public final double currentValue(Atom a) {
        if(phase != a.parentPhase()) {return 0.0;}  //also handles condition that phase contains no atoms
        double pe = 0.0;
        apiUp.reset(a,Iterator.INTER);
        while(apiUp.hasNext()) {
            AtomPair pair = apiUp.next();
            pe += phase.parentSimulation.getPotential(pair).energy(pair);
//            pe += pair.potential.energy(pair);
        }
        apiDown.reset(a,Iterator.INTER);
        while(apiDown.hasNext()) {
            AtomPair pair = apiDown.next();
            pe += phase.parentSimulation.getPotential(pair).energy(pair);
//            pe += pair.potential.energy(pair);
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
            pe += phase.parentSimulation.getPotential(pair).energy(pair);
//            pe += pair.potential.energy(pair);
        }    
      return pe;
    }

    
}
