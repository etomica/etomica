package simulate;
import java.awt.Rectangle;

/**
 * Meter for measurement of the local-molecule-number density inside of a rectangular volume
 */

public class MeterLocalDensityCube extends simulate.MeterLocalDensity
{
    double xCenter, yCenter, halfWidth, halfHeight;
    
    public MeterLocalDensityCube()
    {
        super();
        setBounds(0,0,300,300);
    }

    /**
     Method to compute the volume of the local region where the density is measured
     Volume is determined converting the size of the component (measured in pixels) to simulation units
     */
    public void computeVolume() {
        Rectangle rectangle = getBounds();
        halfWidth = 0.5*rectangle.width / Phase.TO_PIXELS;
        halfHeight = 0.5*rectangle.height / Phase.TO_PIXELS;
        volume = 4 * halfWidth * halfHeight;
        xCenter = rectangle.x/Phase.TO_PIXELS + halfWidth;
        yCenter = rectangle.y/Phase.TO_PIXELS + halfHeight;
    }
    
    /**
     Method that specifies if a molecule is inside the local region where the density is measured
     */
    public boolean contains(Molecule m) {
        double[] r = m.COM();  //molecule center-of-mass
        if(Math.abs(r[0]-xCenter) > halfWidth) {return false;}
        else if(Math.abs(r[1]-yCenter) > halfHeight) {return false;}
        else {return true;}
    }
}
