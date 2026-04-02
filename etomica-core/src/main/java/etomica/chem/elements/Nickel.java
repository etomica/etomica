package etomica.chem.elements;

public class Nickel extends ElementChemical {

    protected Nickel(String symbol) {
        this(symbol, 56.69);
    }

    protected Nickel(String symbol, double mass) {
        super(symbol, mass, 28);
    }

    public static final Nickel INSTANCE = new Nickel("Ni");

    private static final long serialVersionUID = 1L;
}
