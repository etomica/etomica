package etomica;

import etomica.units.Dimension;

/**
 * Meter for the total kinetic energy in a phase
 * Computes total KE by summing values of KE returned by every atom in the phase
 */
 
 /* History of changes
  * 7/03/02  Added non-registering constructor (space argument)
  */
public class MeterKineticEnergy extends MeterScalar
{
    AtomIterator atomIterator;
    
    public MeterKineticEnergy() {
        this(Simulation.instance);
    }
    public MeterKineticEnergy(SimulationElement parent) {
        super(parent);
        setLabel("Kinetic Energy");
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total kinetic energy of molecular motion in a phase");
        return info;
    }


    public Dimension getDimension() {return Dimension.ENERGY;}

    /**
     * This meter needs iterators to do its measurements, so this method overrides the no-op method of AbstractMeter 
     * It obtains the necessary iterators from the phase.
     */
	protected void setPhaseIteratorFactory(IteratorFactory factory) {
       // atomIterator = factory.makeAtomIterator();
       atomIterator = phase.makeAtomIterator();
	}
	
	public void setPhase(Phase p) {
	    super.setPhase(p);
	    atomIterator = p.makeAtomIterator();
	}
	
    public double currentValue() {
        double ke = 0.0;
        atomIterator.reset();
        while(atomIterator.hasNext()) {    //consider doing this with an allAtoms call
            Atom atom = atomIterator.next();
            if(atom.type instanceof AtomType.Wall) continue;
            ke += atom.coord.kineticEnergy();
//            ke += atomIterator.next().coord.kineticEnergy();
        }
        return ke;
    }//end of currentValue
    
/*    public static double currentValue(Phase p) {
        double ke = 0.0;
        for(Atom atom=p.firstAtom(); atom!=null; atom=atom.nextAtom()) {
            if(atom.type instanceof AtomType.Wall) continue;
            ke += atom.coord.kineticEnergy();
        }
        return ke;
    }//end of value*/
}