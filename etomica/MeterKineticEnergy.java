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
    AtomIteratorList atomIterator = new AtomIteratorList();
    
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
	
    public double getDataAsScalar(Phase phase) {
        double ke = 0.0;
        atomIterator.setList(phase.speciesMaster.atomList);
        atomIterator.reset();
        while(atomIterator.hasNext()) {    //consider doing this with an allAtoms call
            Atom atom = atomIterator.next();
            if(atom.type instanceof AtomType.Wall) continue;
            ke += atom.coord.kineticEnergy();
//            ke += atomIterator.next().coord.kineticEnergy();
        }
        return ke;
    }//end of getDataAsScalar
 }