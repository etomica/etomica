package etomica.chem.elements;

public class Oxygen_2 extends ElementChemical {

    protected Oxygen_2(String symbol) {
        this(symbol, 15.9994);
    }

    protected Oxygen_2(String symbol, double mass) {
        super(symbol, mass, 8);
    }

    private static final long serialVersionUID = 1L;
    public static final Oxygen_2 INSTANCE = new Oxygen_2("O_2");
}

