package etomica;
import etomica.units.Dimension;

/**
 * Meter for evaluation of the potential energy in a phase
 * Includes several related methods for computing the potential energy of a single
 * atom or molecule with all neighboring atoms
 *
 * @author David Kofke
 */
 
public class MeterPotentialEnergy extends MeterScalar implements EtomicaElement {
    
    public MeterPotentialEnergy() {
        this(Simulation.instance);
    }
    public MeterPotentialEnergy(Simulation sim) {
        super(sim);
        setLabel("Potential Energy");
        iteratorDirective.includeLrc = true;
        potential = sim.hamiltonian.potential;
        energy = sim.energySum(this);
    }
      
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total intermolecular potential energy in a phase");
        return info;
    }

    public Dimension getDimension() {return Dimension.ENERGY;}
    
    public void setIncludeLrc(boolean b) {
    	iteratorDirective.includeLrc = b;
    }
    public boolean isIncludeLrc() {
    	return iteratorDirective.includeLrc;
    }

    public void setTarget(Atom atom) {
    	singleAtom[0] = atom;
    	iteratorDirective.setTargetAtoms(singleAtom);
    }
    
    public void setTarget(Atom[] atoms) {
    	iteratorDirective.setTargetAtoms(atoms);
    }
 /**
  * Computes total potential energy for phase.
  * Currently, does not include long-range correction to truncation of energy
  */
    public double getDataAsScalar(Phase phase) {
    	energy.zeroSum();
        potential.calculate(phase, iteratorDirective, energy);
        return energy.getSum();
    }
    
    private final IteratorDirective iteratorDirective = new IteratorDirective();
    private final PotentialCalculationEnergySum energy;
    private final PotentialMaster potential;
    private final Atom[] singleAtom = new Atom[1];
    
}//end of MeterPotentialEnergy
