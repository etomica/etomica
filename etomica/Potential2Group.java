package etomica;


public class Potential2Group extends Potential2 implements PotentialGroup {
    
    private final IteratorDirective localDirective = new IteratorDirective();
    private final PotentialCalculation.EnergySum energy = new PotentialCalculation.EnergySum();
    private PotentialLinker first;
    
    public Potential2Group(PotentialGroup parent) {
        super(parent);
    }
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        iterator.reset(id);
 //       localDirective.copy(id);
        while(iterator.hasNext()) {
            AtomPair pair = iterator.next();
   //         localDirective.setBasis(pair.atom1(), pair.atom2());
            for(PotentialLinker link=first; link!=null; link=link.next) {
                if(id.excludes(link.potential)) continue; //see if potential is ok with iterator directive
                link.potential.set(pair.atom1, pair.atom2).calculate(id, pc);
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

    public double energy(AtomPair pair) {
        return calculate(localDirective.set().set(IteratorDirective.BOTH), energy).sum();
    }
        

}//end Potential2Group
    