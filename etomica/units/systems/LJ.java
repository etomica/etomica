package etomica.units.systems;
import etomica.Default;
import etomica.units.BaseUnit;
import etomica.units.Count;
import etomica.units.Decimal;
import etomica.units.Prefix;
import etomica.units.Radian;
import etomica.units.Unit;
import etomica.units.UnitSystem;

/**
 * Lennard-Jones system of units, such that all quantities are made
 * dimensionless with respect to a characteristic size (sigma), energy (epsilon)
 * and mass. Values of sigma, epsilon and mass may be set to any value using the
 * static accessor methods (set/get).  Changes are immediately propagated to the
 * sim- >LJ unit conversion factors.
 */
 
 /* History
  * 03/12/04 (DAK) removed most static features, allowing for different
  * instances
  */
  
public class LJ extends UnitSystem implements java.io.Serializable {
    
    private double sigma = Default.ATOM_SIZE;
    private double epsilon = Default.POTENTIAL_WELL;
    private double mass = Default.ATOM_MASS;
    
    private Mass massUnit = new Mass(sigma, epsilon, mass);
    private Length lengthUnit = new Length(sigma, epsilon, mass);
    private Time timeUnit = new Time(sigma, epsilon, mass);
    private Charge chargeUnit = new Charge(sigma, epsilon, mass);
    private Dipole dipoleUnit = new Dipole(sigma, epsilon, mass);
    private Energy energyUnit = new Energy(sigma, epsilon, mass);
    private Temperature temperatureUnit = new Temperature(sigma, epsilon, mass);
    private Pressure2D pressure2DUnit = new Pressure2D(sigma, epsilon, mass);
    private Pressure pressureUnit = new Pressure(sigma, epsilon, mass);
    private Volume2D volume2DUnit = new Volume2D(sigma, epsilon, mass);
    private Volume volumeUnit = new Volume(sigma, epsilon, mass);
        
    public double getSigma() {return sigma;}
    public void setSigma(double s) {sigma = s; update();}
    public double getEpsilon() {return epsilon;}
    public void setEpsilon(double e) {epsilon = e; update();}
    public double getMass() {return mass;}
    public void setMass(double m) {mass = m; update();}
    
	public Unit quantity() {return Count.UNIT;}
    public Unit fraction() {return Decimal.UNIT;}
	public Unit mass() {return massUnit;}
	public Unit length() {return lengthUnit;}
	public Unit time() {return timeUnit;}
	public Unit angle() {return Radian.UNIT;}
	public Unit charge() {return chargeUnit;}  
	public Unit dipole() {return dipoleUnit;}  
	public Unit energy() {return energyUnit;}
	public Unit temperature() {return energyUnit;}
	public Unit pressure(int D) {return (D==2) ? (Unit)pressure2DUnit : (Unit)pressureUnit;}  
	public Unit volume(int D) {return (D==2) ? (Unit)volume2DUnit : (Unit)volumeUnit;} 

        
    private void update() {
        massUnit.updateConversions(sigma, epsilon, mass);
        lengthUnit.updateConversions(sigma, epsilon, mass);
        timeUnit.updateConversions(sigma, epsilon, mass);
        chargeUnit.updateConversions(sigma, epsilon, mass);
        dipoleUnit.updateConversions(sigma, epsilon, mass);
        energyUnit.updateConversions(sigma, epsilon, mass);
        temperatureUnit.updateConversions(sigma, epsilon, mass);
        pressureUnit.updateConversions(sigma, epsilon, mass);
        pressure2DUnit.updateConversions(sigma, epsilon, mass);
        volumeUnit.updateConversions(sigma, epsilon, mass);
        volume2DUnit.updateConversions(sigma, epsilon, mass);
    }
    
    //  \u03B5 is unicode for epsilon
    //  \u03C3 is unicode for sigma
    //  \u00BD is unicode for "1/2"
    
    public static final class Mass extends BaseUnit.Mass {
        private Mass(double sigma, double epsilon, double mass) {
        	super(1.0, "LJ units", "m", Prefix.NOT_ALLOWED);
        	updateConversions(sigma, epsilon, mass);
        }
        private void updateConversions(double sigma, double epsilon, double mass) {
            setToSimConversionFactor(mass); //conversion from LJ unit to simulation unit
        }
    }

    public static final class Length extends BaseUnit.Length {
        private Length(double sigma, double epsilon, double mass) {
			super(1.0, "LJ units", "\u03C3", Prefix.NOT_ALLOWED);
			updateConversions(sigma, epsilon, mass);
        }
        private void updateConversions(double sigma, double epsilon, double mass) {
			setToSimConversionFactor(sigma); //conversion from LJ unit to simulation unit
        }
    }

    public static final class Time extends BaseUnit.Time {
        private Time(double sigma, double epsilon, double mass) {
			super(1.0, "LJ units", "(\u03B5m\u03C3^2)^\u00BD", Prefix.NOT_ALLOWED);
            updateConversions(sigma, epsilon, mass);
        }
        private void updateConversions(double sigma, double epsilon, double mass) {
			setToSimConversionFactor(1.0/Math.sqrt(epsilon * mass * sigma*sigma)); //conversion from LJ unit to simulation unit
        }
    }

    public static final class Charge extends BaseUnit.Charge {
        private Charge(double sigma, double epsilon, double mass) {
			super(1.0, "LJ units", "(\u03B5\u03C3)^(-\u00BD)", Prefix.NOT_ALLOWED);
            updateConversions(sigma, epsilon, mass);
        }
        private void updateConversions(double sigma, double epsilon, double mass) {
			setToSimConversionFactor(Math.sqrt(epsilon * sigma)); //conversion from LJ unit to simulation unit
        }
    }

    public static final class Dipole extends BaseUnit.Dipole {
        private Dipole(double sigma, double epsilon, double mass) {
			super(1.0, "LJ units", "(\u03B5\u03C3^3)^(-\u00BD)", Prefix.NOT_ALLOWED);
            updateConversions(sigma, epsilon, mass);
        }
        private void updateConversions(double sigma, double epsilon, double mass) {
			setToSimConversionFactor(Math.sqrt(epsilon * sigma*sigma*sigma)); //conversion from LJ unit to simulation unit
        }
    }

    public static final class Energy extends BaseUnit.Energy {
        private Energy(double sigma, double epsilon, double mass) {
			super(1.0, "LJ units", "\u03B5", Prefix.NOT_ALLOWED);
            updateConversions(sigma, epsilon, mass);
        }
        private void updateConversions(double sigma, double epsilon, double mass) {
			setToSimConversionFactor(epsilon); //conversion from LJ unit to simulation unit
        }
    }


    public static final class Temperature extends BaseUnit.Temperature {
        private Temperature(double sigma, double epsilon, double mass) {
			super(1.0, "LJ units", "\u03B5", Prefix.NOT_ALLOWED);
            updateConversions(sigma, epsilon, mass);
        }
        private void updateConversions(double sigma, double epsilon, double mass) {
			setToSimConversionFactor(epsilon); //conversion from LJ unit to simulation unit
        }
    }
    public static final class Pressure extends BaseUnit.Pressure {
        private Pressure(double sigma, double epsilon, double mass) {
			super(1.0, "LJ units", "\u03B5/\u03C3^3", Prefix.NOT_ALLOWED);
            updateConversions(sigma, epsilon, mass);
        }
        private void updateConversions(double sigma, double epsilon, double mass) {
			setToSimConversionFactor(epsilon/(sigma*sigma*sigma)); //conversion from LJ unit to simulation unit
        }
    }

    public static final class Pressure2D extends BaseUnit.Pressure2D {
        private Pressure2D(double sigma, double epsilon, double mass) {
			super(1.0, "LJ units", "\u03B5/\u03C3^2", Prefix.NOT_ALLOWED);
            updateConversions(sigma, epsilon, mass);
        }
        private void updateConversions(double sigma, double epsilon, double mass) {
			setToSimConversionFactor(epsilon/(sigma*sigma)); //conversion from LJ unit to simulation unit
        }
    }

    public static final class Volume extends BaseUnit.Volume {
        private Volume(double sigma, double epsilon, double mass) {
			super(1.0, "LJ units", "\u03C3^3", Prefix.NOT_ALLOWED);
            updateConversions(sigma, epsilon, mass);
        }
        private void updateConversions(double sigma, double epsilon, double mass) {
			setToSimConversionFactor(sigma*sigma*sigma); //conversion from LJ unit to simulation unit
        }
    }

    public static final class Volume2D extends BaseUnit.Volume2D {
        private Volume2D(double sigma, double epsilon, double mass) {
			super(1.0, "LJ units", "\u03C3^2", Prefix.NOT_ALLOWED);
            updateConversions(sigma, epsilon, mass);
        }
        private void updateConversions(double sigma, double epsilon, double mass) {
			setToSimConversionFactor(sigma*sigma); //conversion from LJ unit to simulation unit
        }
    }
}