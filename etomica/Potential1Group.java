package etomica;


public class Potential1Group extends Potential1 implements PotentialGroup {
    
    private final IteratorDirective localDirective = new IteratorDirective();
    private final PotentialCalculation.EnergySum energy = new PotentialCalculation.EnergySum();
    protected PotentialLinker first;
    
    public Potential1Group(PotentialGroup parent) {
        super(parent);
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

    public void addPotential(Potential potential) {
        first = new PotentialLinker(potential, first);
    }

    public double energy(Atom atom) {
        calculate(localDirective.set(atom).set(IteratorDirective.BOTH), energy);
        return energy.sum();
    }

}//end Potential1Group
    