package etomica.chem.elements;

public class Boron extends ElementChemical {

    protected Boron(String symbol) {
        this(symbol, 10.81);
    }

    protected Boron(String symbol, double mass) {
        super(symbol, mass, 5);
    }

    public static final Boron INSTANCE = new Boron("B");

    private static final long serialVersionUID = 1L;
}
