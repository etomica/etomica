package etomica.chem.elements;

public class Neodymium extends ElementChemical {

    protected Neodymium(String symbol) {
        this(symbol, 144.24);
    }
    protected Neodymium(String symbol, double mass) {
        super(symbol, mass, 60);
    }

    public static final Neodymium INSTANCE = new Neodymium("Nd");

    private static final long serialVersionUID = 1L;
}
