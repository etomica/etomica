package simulate;

public abstract class AtomP extends Atom {
    
    public LatticeSite site;
    
    public AtomP(Molecule parent, int index) {
        super(parent, index);
        if(parent != null) {    //null parent indicates atom is used only to generate other atoms
            setStationary(false);
        }
    }
}