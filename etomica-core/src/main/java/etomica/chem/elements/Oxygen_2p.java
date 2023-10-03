package etomica.chem.elements;

public class Oxygen_2p extends ElementChemical {

    protected Oxygen_2p(String symbol) {
        this(symbol, 15.9994);
    }

    protected Oxygen_2p(String symbol, double mass) {
        super(symbol, mass, 8);
    }

    private static final long serialVersionUID = 1L;
    public static final Oxygen_2p INSTANCE = new Oxygen_2p("O_2p");
}

