package etomica.chem.elements;

public class Carbon_Arp extends ElementChemical{
    protected Carbon_Arp(String symbol) {
        this(symbol, 12.0107);
    }

    protected Carbon_Arp(String symbol, double mass) {
        super(symbol, mass, 6);
    }

    public static final Carbon_Arp INSTANCE = new Carbon_Arp("C_Arp");

    private static final long serialVersionUID = 1L;
}
