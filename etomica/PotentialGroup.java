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
    
    public double energy() {
        double sum = 0.0;
        for(PotentialAbstract p=first; p!=null; p=p.next) {
           sum += p.energy();
        }
        return sum;
    }
/*    public double energy(AtomGroup group) {
        double sum = 0.0;
        for(PotentialAbstract p=first; p!=null; p=p.next) {
           sum += p.energy(group);
        }
        return sum;
    }
    */
    public double energy(Atom atom) {
        double sum = 0.0;
        for(PotentialAbstract p=first; p!=null; p=p.next) {
           sum += p.energy(atom);
        }
        return sum;
    }
    
/*    public double energy(AtomGroup[] group) {
        double sum = 0.0;
        iterator.reset(group);
        while(iterator.hasNext()) {
            AtomGroup atomGroup = iterator.next();
            for(PotentialGroup p=first; p!=null; p=p.next) {
                sum += p.energy(atomGroup);
            }
        }
        return sum;
    }
 */           
    public void addPotential(PotentialAbstract potential) {
        if(first == null) first = potential;
        if(last != null) last.next = potential;
 //       potential.previous = last;
        last = potential;
    }
        
}//end of PotentialGroup
    