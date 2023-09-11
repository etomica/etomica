package etomica.chem.elements;

public class Vanadium extends ElementChemical {

	protected Vanadium(String symbol) {this(symbol, 50.94);}

    protected Vanadium(String symbol, double mass) {
        super(symbol, mass, 23);
    }

    public static final Vanadium INSTANCE = new Vanadium("V");

    private static final long serialVersionUID = 1L;
}
