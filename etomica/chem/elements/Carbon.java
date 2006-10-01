package etomica.chem.elements;

/**
 * Class for the Carbon element. 
 *
 * @author Andrew
 */
public class Carbon extends ElementChemical {
    
    protected Carbon(String symbol) {
        this(symbol, 12.0107);
    }
    
    protected Carbon(String symbol, double mass) {
        super(symbol, mass, 6);
    }

    public static final Carbon INSTANCE = new Carbon("C");
}
