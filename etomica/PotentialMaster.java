package etomica;

/**
 * Master potential that sits that the top of the hierarchy of
 * potentials in a simulation.  
 *
 * @author David Kofke
 */
public final class PotentialMaster extends Potential1Group /*implements java.io.Serializable */{
    
    private PotentialLinker first;
    private SpeciesMaster speciesMaster;

    public PotentialMaster(Simulation sim) {
        super(sim);
    }
    
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(id.excludes(link.potential)) continue; //see if potential is ok with iterator directive
            link.potential.calculate(id, pc);
        }//end for
    }//end calculate
        
/*    public final PotentialCalculation.Sum calculate(IteratorDirective id, PotentialCalculation.Sum pa) {
        this.calculate(id, (PotentialCalculation)pa);
        return pa;
    }
                
    public void addPotential(Potential potential) {
        first = new PotentialLinker(potential, first);
        potential.set(speciesMaster);
    }
*/
    //Sets the basis for iteration
    public Potential set(Atom a) {
        speciesMaster = (SpeciesMaster)a;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            link.potential.set(speciesMaster);
        }//end for
        return this;
    }
//    public Potential set(Atom a1, Atom a2) {return null;} //exception?
    public Potential set(Phase p) {return set(p.speciesMaster);}

}//end of PotentialMaster
    