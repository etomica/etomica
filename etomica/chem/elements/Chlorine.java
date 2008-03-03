package etomica.chem.elements;

/**
 * Class for the Chlorine element. 
 *
 * @author Andrew
 */
public class Chlorine extends ElementChemical {

	protected Chlorine(String symbol) {
        this(symbol, 35.453);
    }
    
    protected Chlorine(String symbol, double mass) {
        super(symbol, mass, 17);
    }

    public static final Chlorine INSTANCE = new Chlorine("Cl");
    
	private static final long serialVersionUID = 1L;
}
