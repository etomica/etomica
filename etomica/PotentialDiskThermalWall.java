package simulate;

public class PotentialDiskThermalWall extends PotentialHardDiskWall {
    
    int heatTransfer;

    PotentialDiskThermalWall() {
        this(0.1);
    }

    PotentialDiskThermalWall(double d) {
        super(d);
    }

/**
 * Handler of collision between disk and wall.  Assigns new velocity to disk based on
 * temperature of wall (velocity is explicitly given in terms of temperature, and is 
 * not sampled from the Boltzmann distribution).  Computes change in kinetic energy
 * of particle resulting from collision, and adds this to qAccumulator in wall.
 * Direction of disk after collision is not randomized; only the magnitude of the velocity is rescaled.
 */
    public void bump(AtomHard atom1, AtomHard atom2)  //this needs updating to check for isStationary
    {
        AtomHardDisk disk;
        AtomHardWall wall;
        if(atom2 instanceof AtomHardWall) {
           disk = (AtomHardDisk)atom1;
           wall = (AtomHardWall)atom2;
        }
        else {
           disk = (AtomHardDisk)atom2;
           wall = (AtomHardWall)atom1;
        }
        
        double pSquaredPrior = Space.v1S(disk.p);
        double kPrior = 0.5*pSquaredPrior*disk.rm;
        double pSquaredAfter = disk.mass*wall.getTemperature()*(double)Space.D/Constants.KE2T;  //need to divide by sqrt(m) to get velocity
        double kAfter = 0.5*pSquaredAfter*disk.rm;
        wall.accumulateQ(kAfter-kPrior);
        double pScale = Math.sqrt(pSquaredAfter/pSquaredPrior);
        
        if(wall.isVertical()) {
            disk.p[1] *= pScale;
            disk.p[0] *= -pScale;
        }
        else if(wall.isHorizontal()) {
            disk.p[1] *= -pScale;
            disk.p[0] *= pScale;
        }
    }
}