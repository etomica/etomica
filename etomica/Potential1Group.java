package etomica;


public class Potential1Group extends Potential1 implements PotentialGroup {
    
    private final IteratorDirective localDirective = new IteratorDirective();
    private final PotentialCalculationEnergySum energy = null;//new PotentialCalculation.EnergySum();
    protected PotentialLinker first;
    
    public Potential1Group() {
        this(Simulation.instance.hamiltonian.potential);
    }
    public Potential1Group(PotentialGroup parent) {
        super(parent);
        //overwrite iterator with one that doesn't force looping over leaf atoms
        iterator = new AtomIteratorSequential(false);
    }
    
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        iterator.reset(id);
        localDirective.copy(id);
        while(iterator.hasNext()) {
            Atom atom = iterator.next();
            if(atom == id.atom1()) localDirective.set();
            for(PotentialLinker link=first; link!=null; link=link.next) {
                if(id.excludes(link.potential)) continue; //see if potential is ok with iterator directive
                link.potential.set(atom).calculate(localDirective, pc);
            }//end for
        }//end while
    }//end calculate

    public void addPotential(Potential potential) {
        first = new PotentialLinker(potential, first);
    }

    public double energy(Atom atom) {
        System.out.println("energy method not implemented in Potential1Group");
        System.exit(1);
  //      calculate(localDirective.set(atom).set(IteratorDirective.BOTH), energy);
        return energy.sum();
    }

}//end Potential1Group
    