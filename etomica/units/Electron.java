package simulate.units;
import simulate.Constants;

/**
 * Unit of charge equal to the magnitude of the charge on an electron
 */
public final class Electron extends BaseUnit.Charge {

  /**
   * Convenience instance of this unit to permit unit to be assigned
   * without creating a new instance of it each time
   */
    public static final Electron UNIT = new Electron();
    
    public Electron() {
        to = 4.803e-10*Math.sqrt(Constants.AVOGADRO*1e24*1e-24); //372.7; conversion from (electron/esu)*(g-cm^3/s^2)^(1/2) to (amu-A^3/ps^2)^(1/2)
        from = 1.0/to;
        name = "electron-charge units";
        symbol = "e";   
    }
}