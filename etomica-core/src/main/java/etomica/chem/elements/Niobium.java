package etomica.chem.elements;

public class Niobium extends ElementChemical {

    protected Niobium(String symbol) {
        this(symbol, 92.906);
    }

    protected Niobium(String symbol, double mass) {
        super(symbol, mass, 41);
    }

    public static final Niobium INSTANCE = new Niobium("Nb");

    private static final long serialVersionUID = 1L;
}