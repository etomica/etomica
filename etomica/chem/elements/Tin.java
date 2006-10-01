package etomica.chem.elements;

/**
 * Class for the Tin element.
 *
 * @author Andrew
 */
public class Tin extends ElementChemical {
    
    protected Tin(String symbol) {
        this(symbol, 118.710);
    }
    
    public Tin(String symbol, double mass) {
        super(symbol, mass, 50);
    }

    public static final Tin INSTANCE = new Tin("Sn");
}
