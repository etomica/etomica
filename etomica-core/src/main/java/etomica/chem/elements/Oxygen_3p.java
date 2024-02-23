package etomica.chem.elements;

public class Oxygen_3p extends ElementChemical {

    protected Oxygen_3p(String symbol) {
        this(symbol, 15.9994);
    }

    protected Oxygen_3p(String symbol, double mass) {
        super(symbol, mass, 8);
    }

    private static final long serialVersionUID = 1L;
    public static final Oxygen_3p INSTANCE = new Oxygen_3p("O_3p");
}

