package etomica.chem.elements;

public class Lithium extends ElementChemical {

    protected Lithium(String symbol) {
        this(symbol, 6.94);
    }

    protected Lithium(String symbol, double mass) {
        super(symbol, mass, 3);
    }

    public static final Lithium INSTANCE = new Lithium("Li");

    private static final long serialVersionUID = 1L;
}

