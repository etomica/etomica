package etomica.chem.elements;

public class Xenon  extends ElementChemical{
    protected Xenon(String symbol) {this(symbol, 131.293);}
    protected Xenon(String symbol, double mass) {
        super(symbol, mass, 54);
    }

    public static final Xenon INSTANCE = new Xenon("Xe");

    private static final long serialVersionUID = 1L;
}
