package etomica.chem.elements;

/**
 * Class for the Sodium element. 
 *
 * @author Andrew
 */
public class Sodium extends ElementChemical {

	protected Sodium(String symbol) {
        this(symbol, 22.98976928);
    }
    
    protected Sodium(String symbol, double mass) {
        super(symbol, mass, 11);
    }

    public static final Sodium INSTANCE = new Sodium("Na");
    
	private static final long serialVersionUID = 1L;
}
