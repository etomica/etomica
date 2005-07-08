package etomica;
import etomica.units.Dimension;

/**
 * Parameters are objects associated with an instance of an AtomType class, and
 * are used to hold parameters characteristic of all atoms having the AtomType
 * instance.
 * 
 * @author David Kofke
 */
 
 public abstract class Parameter implements java.io.Serializable {
    
    public static String VERSION = "Parameter:01.11.23";
    
    public interface Size {
        public double getSigma();
        public void setSigma(double s);
        public Dimension getSigmaDimension();
        public Dimension DIMENSION = Dimension.LENGTH;
    }
    
    public interface Energy {
        public double getEpsilon();
        public void setEpsilon(double e);
        public Dimension getEpsilonDimension();
        public Dimension DIMENSION = Dimension.ENERGY;
    }
    
    public interface Mass {
        public double getMass();
        public void setMass(double m);
        public Dimension getMassDimension();
        public Dimension DIMENSION = Dimension.MASS;
    }
    
   
    public interface Source {
        public Parameter makeParameter();
    }
    
 }
