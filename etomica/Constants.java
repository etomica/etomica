package simulate;

public class Constants extends Object {
    
    private Constants() {}   // can't instantiate class
    
    /* Assumed units for I/O
             time: picoseconds
           length: Angstroms
      temperature: Kelvins
             mass: amu
         pressure: bar
             
       Units for internal calculations
            time:  
    */
    
    public static final double KE2T = 10.0/8.314;   //converts kinetic energy (mass*velocity^2/kB) to temperature
    public static final double PV2T = 1e5*1e-30/8.314*6.022e23;  //converts P*V/kB to temperature in Kelvins
}
    
    
    