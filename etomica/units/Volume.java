package etomica.units;

import java.io.ObjectStreamException;

/**
 * Base for all volume units. Simulation unit of volume is A^3
 */
public final class Volume extends Dimension {
    
    public static final Dimension DIMENSION = new Volume();
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1, "cubic Angstroms", "\u00c5^3", Prefix.NOT_ALLOWED);

    private Volume() {
        super("Volume", 3, 0, 0, 0, 0, 0, 0);
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.volume();
    }

   public static Dimension dimension(int D) {
        switch(D) {
            case 3: 
                return Volume.DIMENSION;
            case 2:
                return Area.DIMENSION;
            case 1:
                return Length.DIMENSION;
            default:
                throw new IllegalArgumentException("Volume dimension defined only for D = 1, 2, 3; you gave D = "+D);
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