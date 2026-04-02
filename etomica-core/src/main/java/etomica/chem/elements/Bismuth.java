package etomica.chem.elements;

public class Bismuth extends ElementChemical {
    protected Bismuth(String symbol) {
        this(symbol, 208.98);
    }
    protected Bismuth(String symbol, double mass) {
        super(symbol, mass, 83);
    }
    public static final Bismuth INSTANCE = new Bismuth("Bi");
    private static final long serialVersionUID = 1L;
}
