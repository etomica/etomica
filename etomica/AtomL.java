package simulate;

public abstract class AtomL extends Atom {
    
    public LatticeSite site;
    
    public AtomL(Molecule parent, int index) {
        super(parent, index);
        if(parent != null) {    //null parent indicates atom is used only to generate other atoms
            setStationary(false);
        }
    }
}