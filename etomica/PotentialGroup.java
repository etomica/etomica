package etomica;

public class PotentialGroup extends PotentialAbstract {
    
    private PotentialAbstract first;
    private PotentialAbstract last;
//    private PotentialGroup parentGroup;
    
    public PotentialGroup() {
        this(Simulation.instance);
    }
    public PotentialGroup(Simulation sim) {
        super(sim);
    }
    
    public double energy(IteratorDirective id) {
        double sum = 0.0;
        for(PotentialAbstract p=first; p!=null; p=p.next) {
           sum += p.energy(id);
        }
        return sum;
    }
    
    public void addPotential(PotentialAbstract potential) {
        if(first == null) first = potential;
        if(last != null) last.next = potential;
 //       potential.previous = last;
        last = potential;
    }
        
}//end of PotentialGroup
    