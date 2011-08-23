package etomica.chem.elements;

/**
 * Class for the Helium element.
 */
public class Helium extends ElementChemical {
	
    protected Helium(String symbol) {
        this(symbol, 4.002602);
    }
	
    protected Helium(String symbol, double mass) {
        super(symbol, mass, 1);
    }
    
    private static final long serialVersionUID = 1L;
    public static final Helium INSTANCE = new Helium("He");
}
