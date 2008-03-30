package etomica.chem.elements;

/**
 * Class for the Hydrogen element. 
 *
 * @author zhaofang
 */
public class Hydrogen extends ElementChemical {
	
    protected Hydrogen(String symbol) {
        this(symbol, 1.00794);
    }
	
    protected Hydrogen(String symbol, double mass) {
        super(symbol, mass, 1);
    }
    
    private static final long serialVersionUID = 1L;
    public static final Hydrogen INSTANCE = new Hydrogen("H");
}
