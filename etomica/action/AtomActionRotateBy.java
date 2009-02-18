package etomica.action;

import java.io.Serializable;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.space.ISpace;
import etomica.space3d.RotationTensor3D;

/**
 * 
 * Performs RIGHT-HANDED rotations of an atom about the x-, y-, and z-axes by the specified angles (roll, pitch, and yaw angles, respectively) in radians.
 * 
 * To move all atoms in a molecule (or atom group), wrap an
 * instance of this class in an AtomGroupAction.
 * 
 * @author Katherine Schadel
 */
public class AtomActionRotateBy implements AtomAction, Serializable {
    
    
    
    public AtomActionRotateBy(ISpace space) {

    }
    
    public void actionPerformed(IAtom atom) {
    	   	
        RotationTensor3D rotate = new RotationTensor3D();
        
        // Rotate about the x-axis
        // Use negative of roll because RotationTensor3D is left-handed about x-axis.
        rotate.setAxial(0, -roll);
        rotate.transform(((IAtomPositioned)atom).getPosition());
        
        // Rotate about the y-axis
        // Use negative of pitch because RotationTensor3D is left-handed about x-axis.
        rotate.setAxial(1, -pitch); 
        rotate.transform(((IAtomPositioned)atom).getPosition());
        
        // Rotate about the z-axis
        // Use negative of yaw because RotationTensor3D is left-handed about x-axis.
        rotate.setAxial(2, -yaw); 
        rotate.transform(((IAtomPositioned)atom).getPosition());
           
    }
    
    public void setAngles (double roll, double pitch, double yaw) {
    	this.roll = roll;
    	this.pitch = pitch;
    	this.yaw = yaw;
    }
    
    private static final long serialVersionUID = 1L;
    private double roll;
    private double pitch;
    private double yaw;

}
