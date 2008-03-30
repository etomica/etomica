package etomica.chem.elements;

/**
 * Class for the Sulfur element.
 *
 * @author Andrew
 */
public class Sulfur extends ElementChemical {
    
    protected Sulfur(String symbol) {
        this(symbol, 32.065);
    }
    
    public Sulfur(String symbol, double mass) {
        super(symbol, mass, 16);
    }

    private static final long serialVersionUID = 1L;
    public static final Sulfur INSTANCE = new Sulfur("S");
}
