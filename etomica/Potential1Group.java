package etomica;


public class Potential1Group extends Potential1 {
    
    private final IteratorDirective localDirective = new IteratorDirective();
    private final PotentialCalculation.EnergySum energy = new PotentialCalculation.EnergySum();
    private PotentialLinker first;
    
    public Potential1Group(Simulation sim) {
        super(sim);
    }
    
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        iterator.reset(id);
//        localDirective.copy(id);
        while(iterator.hasNext()) {
            Atom atom = iterator.next();
//            localDirective.set(atom);
            for(PotentialLinker link=first; link!=null; link=link.next) {
                if(id.excludes(link.potential)) continue; //see if potential is ok with iterator directive
                link.potential.set(atom).calculate(id, pc);
            }//end for
        }//end while
    }//end calculate

    public final PotentialCalculation.Sum calculate(IteratorDirective id, PotentialCalculation.Sum pa) {
        this.calculate(id, (PotentialCalculation)pa);
        return pa;
    }
        
    public void addPotential(Potential potential) {
        first = new PotentialLinker(potential, first);
    }

    public double energy(Atom atom) {
        return calculate(localDirective.set().set(IteratorDirective.BOTH), energy).sum();
    }
        

}//end Potential1Group
    