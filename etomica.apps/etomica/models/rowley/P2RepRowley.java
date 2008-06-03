package etomica.models.rowley;

import etomica.EtomicaInfo;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.ISpace;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Purely repulsive potential from Rowley et al (2006) used for interactions between satellite sites, X.
 * These fudge sites are used to represent a region of high electron density belonging to the oxygen atom of simple alcohols.    
 *
 * K.R. Schadel
 * May 2008
 */
public final class P2RepRowley extends Potential2SoftSpherical {

    public P2RepRowley (ISpace space) {
        this(space, 1.0, 1.0);
    }
    
    public P2RepRowley (ISpace space, double BXX, double CXX) {
        super(space);
        setBXX(BXX);
        setCXX(CXX);

    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Purely repulsive potential from Rowley et al 2006");
        return info;
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
    	
    	double r = Math.sqrt(r2);
    	
    	return BXX*Math.exp(-CXX*r);
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
    	
    	return 0;
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
 	
        return 0;
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
    	
        return 0;  
    }


    public double getBXX() {return BXX;}
 
    public final void setBXX(double eps) {
    	BXX = eps;
    }
    public Dimension getBXXDimension() {return Energy.DIMENSION;}
   
    
    public double getCXX() {return CXX;}

    public final void setCXX(double a) {
    	CXX = a;
    }
    
    private static final long serialVersionUID = 1L;
    private double BXX;
    private double CXX;
 
}

