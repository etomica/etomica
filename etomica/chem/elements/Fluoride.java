package etomica.chem.elements;

/**
 * Class for the Sulfur element.
 *
 * @author Andrew
 */
public class Fluoride extends ElementChemical {
    
    protected Fluoride(String symbol) {
        this(symbol, 18.998404);
    }
    
    public Fluoride(String symbol, double mass) {
        super(symbol, mass, 9);
    }

    private static final long serialVersionUID = 1L;
    public static final Fluoride INSTANCE = new Fluoride("F");
}
