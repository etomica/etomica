package etomica;

 /* History of changes
  * 8/14/02 (DAK) introduced removePotential method and modified addPotential.
  */

public class Potential1Group extends Potential1 implements PotentialGroup {
    
    private final IteratorDirective localDirective = new IteratorDirective();
    private final PotentialCalculationEnergySum energy = null;//new PotentialCalculation.EnergySum();
    protected PotentialLinker first;
    private final Atom[] atoms = new Atom[1];
    
    public Potential1Group() {
        this(Simulation.instance.hamiltonian.potential);
    }
    public Potential1Group(PotentialGroup parent) {
        super(parent);
    }
    
    public boolean contains(Potential potential) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential.equals(potential)) return true;
        }//end for
        return false;
    }
    
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        iterator.reset(id);
        localDirective.copy(id);
        while(iterator.hasNext()) {
            atoms[0] = iterator.next();
            if(atoms[0] == id.atom1()) localDirective.set();
            for(PotentialLinker link=first; link!=null; link=link.next) {
                if(id.excludes(link.potential)) continue; //see if potential is ok with iterator directive
                link.potential.set(atoms).calculate(localDirective, pc);
            }//end for
        }//end while
    }//end calculate

    /**
     * Adds the given potential to this group, but should not be called directly.  Instead,
     * this method is invoked by the setParentPotential method (or more likely, 
     * in the constructor) of the given potential.  
     * Part of the PotentialGroup interface (so must be defined with public access).
     */
    public void addPotential(Potential potential) {
        if(potential == null || potential.parentPotential() != this) 
            throw new IllegalArgumentException("Improper call to addPotential; should use setParentPotential method in child instead of addPotential method in parent group"); 
        first = new PotentialLinker(potential, first);
    }

    /**
     * Removes given potential from the group.  No error is generated if
     * potential is not in group.
     */
    public void removePotential(Potential potential) {
        PotentialLinker previous = null;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential == potential) {//found it
                if(previous == null) first = link.next;  //it's the first one
                else previous.next = link.next;          //it's not the first one
                return;
            }
            previous = link;
        }
    }

    public double energy(Atom atom) {
        throw new RuntimeException("energy method not implemented in Potential1Group");
  //      calculate(localDirective.set(atom).set(IteratorDirective.BOTH), energy);
  //      return energy.sum();
    }

}//end Potential1Group
    