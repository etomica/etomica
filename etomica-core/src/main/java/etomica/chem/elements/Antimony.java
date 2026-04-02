package etomica.chem.elements;

public class Antimony extends ElementChemical {
    protected Antimony(String symbol) {
        this(symbol, 121.76);
    }
    protected Antimony(String symbol, double mass) {
        super(symbol, mass, 51);
    }
    public static final Antimony INSTANCE = new Antimony("Sb");
    private static final long serialVersionUID = 1L;
}
