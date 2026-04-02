package etomica.chem.elements;

public class Gadolinium extends ElementChemical {

    protected Gadolinium(String symbol) {
        this(symbol, 157.25);
    }

    protected Gadolinium(String symbol, double mass) {
        super(symbol, mass, 64);
    }

    public static final Gadolinium INSTANCE = new Gadolinium("Gd");

    private static final long serialVersionUID = 1L;
}
