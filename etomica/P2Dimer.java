package etomica;

/**
 * Dissociable dimer potential
 */
 
public class P2Dimer extends Potential2 {
    
    Potential bonded;
    Potential oneRadical;
    Potential bothRadicals;
    Potential bothSaturated;
    
    public P2Dimer() {
        this(Simulation.instance);
    }
    
    public P2Dimer(Simulation sim) {
        super(sim);
    }
    
    public P2Dimer(Simulation sim, 
        Potential bonded, Potential oneRadical, Potential bothRadicals, Potential bothSaturated) {
            super(sim);
            this.bonded = bonded;
            this.oneRadical = oneRadical;
            this.bothRadicals = bothRadicals;
            this.bothSaturated = bothSaturated;
   }
    
    public Potential getPotential(Atom a1, Atom a2) {
        if (a1.atomLink[0] == null) {
            if(a2.atomLink[0] == null) { 
                return bothRadicals;
            }
            else {
                return oneRadical;
            }
        }
        else { //a1 is saturated
            if(a2.atomLink[0] == null) {
                return oneRadical;
            }
            if (a1.atomLink[0].atom() == a2) return bonded;
            else return bothSaturated;
        }
    }
}//end of class        
    
    