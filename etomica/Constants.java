package simulate;
import java.awt.Color;

public class Constants extends Object {
    
    private Constants() {}   // can't instantiate class
    
    /* Assumed units for I/O
             time: picoseconds (ps)
           length: Angstroms   (A)
      temperature: Kelvins     (K)
             mass: amu
         pressure: bar
             
       Units for internal calculations (simulation units)
            time: ps
          length: Angstroms
            mass: amu
    */
    
    // Scaling factor from simulation length scale to Angstroms (for temporary use) (ha ha!)
    public static final double SCALE = 35.;
    public static final double DEPTH = 3.5;
    
    public static final double AVOGADRO = 6.022e23;
    
    public static final double KE2T = 10.0/8.314;   //converts kinetic energy (mass*velocity^2/kB) to temperature
    public static final double PV2T = 1e5*1e-30/8.314*AVOGADRO;  //converts P*V/kB to temperature in Kelvins
    public static final double BAR2SIM = 1e5*1000*AVOGADRO/1e10/1e24;  //converts bar to amu/A-ps^2 (pressure)
    public static final double MOL_PER_LITER2SIM = 1e3*Constants.AVOGADRO/1e30; //converts mol/liter to molecules/A^3
    
    
    public static double G = 9.8*1e10/1e24;  //acceleration of gravity (on Earth), in A/ps^2

  /* Colors adopted in the web textbook on molecular simulation */
    public static final Color KHAKI = new Color(153,153,102);
    public static final Color DARK_KHAKI = new Color(102,102,051);
    public static final Color BRIGHT_RED = new Color(153,000,000);
    public static final Color DARK_RED = new Color(102,000,000);
    public static final Color BLUSH = new Color(153,102,102);
    public static final Color TAN = new Color(204,204,153);
}
    
    
    