package etomica.chem.elements;

public class Palladium extends ElementChemical {

	protected Palladium(String symbol) {this(symbol, 106.4);}

    protected Palladium(String symbol, double mass) {
        super(symbol, mass, 46);
    }

    public static final Palladium INSTANCE = new Palladium("Pd");

    private static final long serialVersionUID = 1L;
}
