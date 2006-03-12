package etomica.units;

/**
 * Unit converting between simulation unit of length and pixels rendered in an
 * image.  Conversion factor is set a construction.
 */
public class Pixel extends SimpleUnit {
	
    /**
     * Constructs with default of 10 pixels per Angstrom.
     */
	public Pixel() {
		this(10.0);
	}
	/**
	 * Constructor for Pixel.  Argument is factor to convert from simulation
	 * units to pixels.  Default is 10.
	 */
	public Pixel(double toPixels) {
		super(Length.DIMENSION, 
			1.0/toPixels, //need to pass superclass conversion from pixels to simulation units
			"pixel","px", Prefix.NOT_ALLOWED
			);
        this.toPixels = toPixels;
	}
    
    public double toPixels() {
        return toPixels;
    }
    
    /**
     * Required because Units are supposed to be singletons
     */
    private Object readResolve() {
        return new Pixel(this.toPixels());
    }

    private final double toPixels;
        
    private static final long serialVersionUID = 1;

}
