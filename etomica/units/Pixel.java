package etomica.units;

import etomica.units.BaseUnit.Length;

/**
 * Unit converting between simulation unit of length and pixels rendered in an
 * image.  Default value is 10 pixels per Angstrom.
 */
public class Pixel extends Length {
	
	/**
	 * Class used for common definition of pixel unit, but not a singleton.
	 */
	public static Pixel UNIT = new Pixel();

	public Pixel() {
		this(10.0);
	}
	/**
	 * Constructor for Pixel.  Argument is factor to convert from simulation
	 * units to pixels.  Default is 10.
	 */
	public Pixel(double toPixels) {
		super(
			1.0/toPixels, //need to pass superclass conversion from pixels to simulation units
			"Pixel", 
			"px", 
			Prefix.NOT_ALLOWED
			);
	}
}
