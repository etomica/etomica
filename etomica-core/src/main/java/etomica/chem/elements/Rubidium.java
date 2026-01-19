package etomica.chem.elements;

public class Rubidium extends ElementChemical {

    protected Rubidium(String symbol) {
        this(symbol, 85.4678);
    }

    protected Rubidium(String symbol, double mass) {
        super(symbol, mass, 37);
    }

    public static final Rubidium INSTANCE = new Rubidium("Rb");

    private static final long serialVersionUID = 1L;
}