package etomica.paracetamol;

import etomica.chem.elements.ElementChemical;

/**
 * Class for the Hydrogen element. 
 *
 * Hp is just designated to be used for generating the CONFIG file for
 * the DL_MULTI package
 *
 */
public class HydrogenP extends ElementChemical {


	protected HydrogenP(String symbol) {
        this(symbol, 1.00794);
    }
	
    protected HydrogenP(String symbol, double mass) {
        super(symbol, mass, 1);
    }
    
    public static final HydrogenP INSTANCE = new HydrogenP("HP");
    private static final long serialVersionUID = 1L;
}
