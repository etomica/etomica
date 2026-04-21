package etomica.chem.elements;

public class Cobalt extends ElementChemical {

	protected Cobalt(String symbol) {this(symbol, 58.93);}

    protected Cobalt(String symbol, double mass) {
        super(symbol, mass, 27);
    }

    public static final Cobalt INSTANCE = new Cobalt("Co");

    private static final long serialVersionUID = 1L;
}
