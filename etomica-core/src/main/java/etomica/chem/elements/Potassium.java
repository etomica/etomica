package etomica.chem.elements;

public class Potassium extends ElementChemical {

    protected Potassium(String symbol) {
        this(symbol, 22.98976928);
    }

    protected Potassium(String symbol, double mass) {
        super(symbol, mass, 11);
    }

    public static final Potassium INSTANCE = new Potassium("K");

    private static final long serialVersionUID = 1L;
}
