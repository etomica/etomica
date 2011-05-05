package etomica.units;
import etomica.util.Constants;

public final class Hartree extends SimpleUnit {

    /**
     * Hartree unit of energy, corresponding to the approximate energy of
     * hydrogen in its ground state.
     */
    public static final Hartree UNIT = new Hartree();
  
    private Hartree() {
        super(Energy.DIMENSION, 
        	2*Constants.PLANCK_H*Constants.RYDBERG_R*Constants.LIGHT_SPEED,  //262549.9617098284
        	"hartrees", "Ha", Prefix.ALLOWED
        	);
    }
    
    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton UNIT
     */
    private Object readResolve() {
        return UNIT;
    }

    private static final long serialVersionUID = 1;

}