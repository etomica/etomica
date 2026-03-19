package etomica.chem.elements;

public class Molybdenum extends ElementChemical {

	protected Molybdenum(String symbol) {this(symbol, 95.96);}

    protected Molybdenum(String symbol, double mass) {
        super(symbol, mass, 42);
    }

    public static final Molybdenum INSTANCE = new Molybdenum("Mo");

    private static final long serialVersionUID = 1L;
}
