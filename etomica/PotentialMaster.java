package etomica;

/**
 * Master potential that oversees all other potentials in the Hamiltonian.  
 *
 * @author David Kofke
 */
public final class PotentialMaster extends Potential1Group {
    
    public String getVersion() {return "PotentialMaster:01.07.23/"+Potential1Group.VERSION;}

    private SpeciesMaster speciesMaster;

    public PotentialMaster(Simulation sim) {
        super(sim);
    }
    
    //should build on this to do more filtering of potentials based on directive
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(id.excludes(link.potential)) continue; //see if potential is ok with iterator directive
            link.potential.calculate(id, pc);
        }//end for
    }//end calculate
        
    public final PotentialCalculation.Sum calculate(IteratorDirective id, PotentialCalculation.Sum pa) {
        this.calculate(id, (PotentialCalculation)pa);
        return pa;
    }
        
    /**
     * Sets the basis for iteration of atoms by all potentials.
     * The given parameter must be an instance of SpeciesMaster.
     */
    public Potential set(Atom a) {
        speciesMaster = (SpeciesMaster)a;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            link.potential.set(speciesMaster);
        }//end for
        return this;
    }

    /**
     * Sets the potentials to iterate on atoms in the given phase.
     */
    public PotentialMaster set(Phase p) {set(p.speciesMaster); return this;}

}//end of PotentialMaster
    