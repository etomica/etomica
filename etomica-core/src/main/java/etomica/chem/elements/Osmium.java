package etomica.chem.elements;

public class Osmium  extends ElementChemical {

    protected Osmium(String symbol) {
        this(symbol, 190.23);
    }

    protected Osmium(String symbol, double mass) {
        super(symbol, mass, 21);
    }

    public static final Osmium INSTANCE = new Osmium("76");

    private static final long serialVersionUID = 1L;
}
