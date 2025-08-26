package etomica.chem.elements;

public class Zirconium extends ElementChemical {

	protected Zirconium(String symbol) {this(symbol, 91.22);}

    protected Zirconium(String symbol, double mass) {
        super(symbol, mass, 40);
    }

    public static final Zirconium INSTANCE = new Zirconium("Zr");

    private static final long serialVersionUID = 1L;
}
