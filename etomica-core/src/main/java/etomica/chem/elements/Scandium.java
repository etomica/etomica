package etomica.chem.elements;

public class Scandium extends ElementChemical {

    protected Scandium(String symbol) {
        this(symbol, 44.955);
    }

    protected Scandium(String symbol, double mass) {
        super(symbol, mass, 21);
    }

    public static final Scandium INSTANCE = new Scandium("Sc");

    private static final long serialVersionUID = 1L;
}
