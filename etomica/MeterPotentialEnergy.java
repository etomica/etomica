package simulate;
import simulate.units.Dimension;

/**
 * Meter for evaluation of the potential energy in a phase
 * Includes several related methods for computing the potential energy of a single
 * atom or molecule with all neighboring atoms
 *
 * @author David Kofke
 */
 
public class MeterPotentialEnergy extends Meter
{
    /**
    * Iterator for a molecule with all other molecules in phase
    */
    AtomPair.Iterator.MP iteratorMP;
    /**
    * Iterator for all pairs in the phase
    */
    AtomPair.Iterator iteratorAllPairs;
    /**
    * Iterator for all atoms in the phase
    */
    Atom.Iterator iteratorAllAtoms;
    /**
    * Iterator for an atom and all neighbor atoms uplist of it
    */
    AtomPair.Iterator apiUp;
    /**
    * Iterator for an atom and all neighbor atoms downlist of it
    */
    AtomPair.Iterator apiDown;
      
    public MeterPotentialEnergy() {
        this(Simulation.instance);
    }
    public MeterPotentialEnergy(Simulation sim) {
        super(sim);
        setLabel("Potential Energy");
    }
      
    /**
     * Declaration that this meter does not use the boundary object of phase when making its measurements
     */
    public final boolean usesPhaseBoundary() {return false;}
    /**
     * Declaration that this meter does the iteratorFactory of phase when making its measurements
     */
    public final boolean usesPhaseIteratorFactory() {return true;}

    public Dimension getDimension() {return Dimension.ENERGY;}
      
    /**
     * This meter needs iterators to do its measurements, so this method overrides the no-op method of AbstractMeter 
     * It obtains the necessary iterators from the phase.
     */
	protected void setPhaseIteratorFactory(IteratorFactory factory) {
        iteratorMP = new AtomPair.Iterator.MP(factory);
        iteratorAllPairs = factory.makeAtomPairIteratorAll();
        iteratorAllAtoms = factory.makeAtomIterator();
        apiUp = factory.makeAtomPairIteratorUp();
        apiDown = factory.makeAtomPairIteratorDown();
	}
	
 /**
  * Computes total potential energy for all atom pairs in phase
  * Returns infinity (MAX_VALUE) as soon as overlap is detected
  * Currently, does not include long-range correction to truncation of energy
  */
    public final double currentValue() {
                        //could make special case for single species, monatomic
        double pe = 0.0;
        //Add in energy of fields with atoms
        for(PotentialField field=phase.firstField(); field!=null; field=field.nextField()) {
            pe += field.energy();  
            if(pe == Double.MAX_VALUE) return Double.MAX_VALUE;
        }

        //Add in all pair interactions
        iteratorAllPairs.reset();
        while(iteratorAllPairs.hasNext()) {
            AtomPair pair = iteratorAllPairs.next();
            double energy = parentSimulation().getPotential(pair).energy(pair);
            if(energy == Double.MAX_VALUE) return Double.MAX_VALUE;
            pe += energy;
        }
//        return pe + LRC();  //put this back when LRC is corrected
        return pe;
    }

    /**
     * Method is currently under development
     * @return full (extensive) long-range correction to the (non-electrostatic) energy
     */
    //This method is untested
    //Gives full (extensive) long-range correction to the energy
    
    //this needs fixing.  Need two separate representative molecules, but now gets two handles to same molecule
    public final double LRC() {
        double lrc = 0.0;
        double V = phase.volume();
        for(Species.Agent s1=phase.firstSpecies(); s1!=null; s1=s1.nextSpecies()) { //loop over species
            int n1 = s1.moleculeCount();    //number of species-1 molecules in phase
            Molecule m1 = s1.parentSpecies().getMolecule();  //a prototype molecule
            m1.atomIterator.reset();
            while(m1.atomIterator.hasNext()) { //loop over all atoms of a molecule of species 1
                Atom a1 = m1.atomIterator.next();
                Molecule m2 = s1.parentSpecies().getMolecule();  //another
                m2.atomIterator.reset();
                while(m2.atomIterator.hasNext()) { //another loop over atoms of species 1 molecule
                    Atom a2 = m2.atomIterator.next();
                    lrc += parentSimulation().getPotential(a1,a2).energyLRC(n1,n1,V);
                }
                for(Species.Agent s2=s1.nextSpecies(); s2!=null; s2=s2.nextSpecies()) { //loop over other species
                    int n2 = s2.moleculeCount();
                    m2 = s2.parentSpecies().getMolecule();  //molecule of the species
                    m2.atomIterator.reset();
                    while(m2.atomIterator.hasNext()) { //loop over atoms of species 2 molecule
                        Atom a2 = m2.atomIterator.next();
                        lrc += 2.0*parentSimulation().getPotential(a1,a2).energyLRC(n1,n2,V);
                    }
                }
            }
        }
        return lrc;
    }

  /**
   * Computes potential energy for atom due to interactions with fields and other atoms in phase.  
   * Returns zero if atom is not in phase.
   * Returns infinity as soon as overlap is encountered
   */
    public final double currentValue(Atom a) {
        if(phase != a.parentPhase()) {return 0.0;}  //also handles condition that phase contains no atoms
        double pe = 0.0;
        //Add in energy of fields with atom; it is responsibility of field to check if this is an affected atom
        for(PotentialField field=phase.firstField(); field!=null; field=field.nextField()) {
            pe += field.energy(a);  
            if(pe == Double.MAX_VALUE) return Double.MAX_VALUE;
        }
        apiUp.reset(a);
        while(apiUp.hasNext()) {
            AtomPair pair = apiUp.next();
            double energy = parentSimulation().getPotential(pair).energy(pair);
            if(energy == Double.MAX_VALUE) return Double.MAX_VALUE;
            pe += energy;
        }
        apiDown.reset(a);
        while(apiDown.hasNext()) {
            AtomPair pair = apiDown.next();
            double energy = parentSimulation().getPotential(pair).energy(pair);
            if(energy == Double.MAX_VALUE) return Double.MAX_VALUE;
            pe += energy;
        }
        return pe;
    }


 /**
  * Computes and returns the total intermolecular energy of a molecule with all molecules in phase
  * Molecule may or may not be in the phase.
  * Does not include intramolecular energy of molecule or interactions with potential-fields
  * Should be checked more carefully regarding handling of inter/intra-molecular atom pairs
  */
    public final double currentValue(Molecule m) {

        iteratorMP.reset(m);  //should take a better look at functioning of iteratorMP
        double pe = 0.0;
        while(iteratorMP.hasNext()) {
            AtomPair pair = iteratorMP.next();
            double energy = parentSimulation().getPotential(pair).energy(pair);
            if(energy == Double.MAX_VALUE) return Double.MAX_VALUE;
            pe += energy;
        }    
        return pe;
    }

    
}
