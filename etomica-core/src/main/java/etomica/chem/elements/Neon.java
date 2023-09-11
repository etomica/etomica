package etomica.chem.elements;

public class Neon extends ElementChemical{
    protected Neon(String symbol) {
        this(symbol, 20.1797);
    }

    protected Neon(String symbol, double mass) {
        super(symbol, mass, 10);
    }

    public static final Neon INSTANCE = new Neon("Ne");

    private static final long serialVersionUID = 1L;
}
