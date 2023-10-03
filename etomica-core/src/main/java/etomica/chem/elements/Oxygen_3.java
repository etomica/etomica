package etomica.chem.elements;

public class Oxygen_3 extends ElementChemical {

    protected Oxygen_3(String symbol) {
        this(symbol, 15.9994);
    }

    protected Oxygen_3(String symbol, double mass) {
        super(symbol, mass, 8);
    }

    private static final long serialVersionUID = 1L;
    public static final Oxygen_3 INSTANCE = new Oxygen_3("O_3");
}

