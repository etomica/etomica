package etomica.chem.elements;

/**
 * Class for the Silver element.
 *
 * @author Andrew
 */
public class Silver extends ElementChemical {
    
    protected Silver(String symbol) {
        this(symbol, 107.8682);
    }
    
    public Silver(String symbol, double mass) {
        super(symbol, mass, 47);
    }

    public static final Silver INSTANCE = new Silver("Ag");
}
