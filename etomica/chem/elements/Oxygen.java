package etomica.chem.elements;

/**
 * Class for the Oxygen element. 
 *
 * @author zhaofang
 */
public class Oxygen extends ElementChemical {
    
    protected Oxygen(String symbol) {
        this(symbol, 15.9994);
    }
    
    protected Oxygen(String symbol, double mass) {
        super(symbol, mass, 8);
    }

    public static final Oxygen INSTANCE = new Oxygen("O");
}
