package etomica.chem.elements;

public class Lanthanum extends ElementChemical {

    protected Lanthanum(String symbol) {
        this(symbol, 138.905);
    }

    protected Lanthanum(String symbol, double mass) {
        super(symbol, mass, 57);
    }

    public static final Lanthanum INSTANCE = new Lanthanum("La");

    private static final long serialVersionUID = 1L;
}
