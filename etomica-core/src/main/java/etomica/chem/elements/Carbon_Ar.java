package etomica.chem.elements;

public class Carbon_Ar extends ElementChemical{
    protected Carbon_Ar(String symbol) {
        this(symbol, 12.0107);
    }

    protected Carbon_Ar(String symbol, double mass) {
        super(symbol, mass, 6);
    }

    public static final Carbon_Ar INSTANCE = new Carbon_Ar("C_Ar");

    private static final long serialVersionUID = 1L;
}
