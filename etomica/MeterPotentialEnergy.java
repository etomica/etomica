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
    
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationEnergySum energy;
    private final PotentialMaster potential;
    
    public MeterPotentialEnergy() {
        this(Simulation.instance);
    }
    public MeterPotentialEnergy(Simulation sim) {
        super(sim);
        setLabel("Potential Energy");
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = true;
        potential = sim.hamiltonian.potential;
        energy = sim.energySum(this);
    }
      
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total intermolecular potential energy in a phase");
        return info;
    }

    public Dimension getDimension() {return Dimension.ENERGY;}
    
    /**
     * Iterator directive specifies the atoms for which the potential energy is measured.
     * Default value measures all atoms in phase.
     */
    public void setIteratorDirective(IteratorDirective directive) {iteratorDirective = directive;}
    
    /**
     * Accessor method for iterator directive.
     */
    public IteratorDirective getIteratorDirective() {return iteratorDirective;}
      
 /**
  * Computes total potential energy for all atom pairs in phase.
  * Currently, does not include long-range correction to truncation of energy
  */
    public final double currentValue() {
  //      int nPairs = (int)potential.set(phase).calculate(iteratorDirective, pairSum.reset()).sum();
  //      System.out.println("pair count: "+nPairs);
  //    return nPairs;
  /*      AtomIterator moleculeIterator = phase.makeMoleculeIterator();
        moleculeIterator.reset();
        potential.set(phase);
        iteratorDirective.set(IteratorDirective.BOTH);
        double sum = 0.0;
        while(moleculeIterator.hasNext()) {
            Atom molecule = moleculeIterator.next();
            sum += potential.calculate(iteratorDirective.set(molecule), energy.reset()).sum();
        }*/
        //System.out.println();
        //System.out.println("  begin in MeterPotentialEnergy");
    //    iteratorDirective.set(IteratorDirective.UP);
        double dbv = potential.calculate(phase, iteratorDirective.set(), energy.reset()).sum();
    //    System.out.println("energies: "+0.5*sum+"  "+dbv);
        //System.out.println("  end in MeterPotentialEnergy");
        //System.out.println();
        return dbv;//potential.set(phase).calculate(iteratorDirective, energy.reset()).sum();
    }
  //  PotentialCalculationPairSum pairSum = new PotentialCalculationPairSum();
    
}//end of MeterPotentialEnergy
