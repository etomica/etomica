package etomica.chem.elements;

public class Arsenic extends ElementChemical {
    protected Arsenic(String symbol) {
        this(symbol, 74.92);
    }
    protected Arsenic(String symbol, double mass) {
        super(symbol, mass, 33);
    }
    public static final Arsenic INSTANCE = new Arsenic("As");
    private static final long serialVersionUID = 1L;
}
