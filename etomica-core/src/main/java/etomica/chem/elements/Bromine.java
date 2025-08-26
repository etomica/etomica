package etomica.chem.elements;

public class Bromine extends ElementChemical {

    protected Bromine(String symbol) {
        this(symbol, 79.90);
    }

    protected Bromine(String symbol, double mass) {
        super(symbol, mass, 35);
    }

    public static final Argon INSTANCE = new Argon("Br");

    private static final long serialVersionUID = 1L;

}

