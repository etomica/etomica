package etomica.chem.elements;

/**
 * Class for the Copper element.
 *
 * @author Andrew
 */
public class Copper extends ElementChemical {
    
    protected Copper(String symbol) {
        this(symbol, 63.546);
    }
    
    public Copper(String symbol, double mass) {
        super(symbol, mass, 29);
    }

    public static final Copper INSTANCE = new Copper("Cu");
}
