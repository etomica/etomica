package etomica.chem.elements;

public class Aluminium extends ElementChemical {
    protected Aluminium(String symbol) {
        this(symbol, 26.98);
    }
    protected Aluminium(String symbol, double mass) {
        super(symbol, mass, 13);
    }
    public static final Aluminium INSTANCE = new Aluminium("Al");
    private static final long serialVersionUID = 1L;
}
