package etomica.chem.elements;

public class Praseodymium extends ElementChemical {

    protected Praseodymium(String symbol) {
        this(symbol, 140.91);
    }

    protected Praseodymium(String symbol, double mass) {
        super(symbol, mass, 59);
    }

    public static final Praseodymium INSTANCE = new Praseodymium("Pr");

    private static final long serialVersionUID = 1L;
}
