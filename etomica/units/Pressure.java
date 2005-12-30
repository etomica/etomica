package etomica.units;

import java.io.ObjectStreamException;

/**
 * Simulation unit of (3D) pressure is (D-A/ps^2)/A^2 = D/(A-ps^2)
 */
public final class Pressure extends Dimension {

    public static final Dimension DIMENSION = new Pressure();
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1, "sim pressure units", "D/(\u00c5-ps^2)", Prefix.NOT_ALLOWED);

    private Pressure() {
        super("Pressure", -1, 1, -2);
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.pressure();
    }

    public static Dimension dimension(int D) {
        switch(D) {
            case 3: 
                return Pressure.DIMENSION;
            case 2:
                return Pressure2D.DIMENSION;
//            case 1:
//                return Length.DIMENSION;
            default:
                throw new IllegalArgumentException("Pressure dimension defined only for D = 2, 3; you gave D = "+D);
        }
    }

    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton DIMENSION
     */
    private Object readResolve() throws ObjectStreamException {
        return DIMENSION;
    }

    private static final long serialVersionUID = 1;

}