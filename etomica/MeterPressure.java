package simulate;

import simulate.*;
import java.beans.*;
import java.awt.*;

// This class needs some serious re-thinking

public class MeterPressure extends simulate.Meter
{
    private P2 momentumP2;
    private SpeciesMeterWalls momentumSpecies;
    private final double vScale = Constants.SCALE*Constants.SCALE*Constants.SCALE;
    private double momentumSum = 0.0;
    private double timeSum = 0.0;
    double diameter = 0.15;

    public MeterPressure()
    {
        super();
        setLabel("Pressure (bar)");
    }
    
    public void initialize() {
        momentumP2 = new P2();
        momentumP2.setSpecies1Index(0);
        momentumP2.setSpecies2Index(1);
        phase.add(momentumP2);
        momentumSpecies = new SpeciesMeterWalls();
        momentumSpecies.setSpeciesIndex(1);
        phase.add(momentumSpecies);
        phase.nMoleculeTotal -= 2;
        phase.nAtomTotal -= 2;
    }

    public void integrationIntervalAction(IntegrationIntervalEvent evt) {
        timeSum += evt.integrator.drawTimeStep * evt.integrator.integrationInterval;
        updateStatistics(phase);
    }

    public double currentValue()
    {
        double flux = 0.5*momentumSum*Constants.SCALE/(timeSum * Constants.SCALE * Constants.DEPTH);
        momentumSum = 0.0;
        timeSum = 0.0;
        return flux/Constants.BAR2SIM;
    }
    
 class P2 extends Potential2 {
    public P2() {
      potential = new Potential[1][1];
      potential[0][0] = new PotentialMomentumFlux();
    }
    public final Potential getPotential(Atom a1, Atom a2) {return potential[0][0];}
  }
    

  class SpeciesMeterWalls extends SpeciesWalls {
    SpeciesMeterWalls() {
        super(2,1);
    }
    public void initializeMolecules() {
      setAngle(90);
      setVisible(false);
    }
    public void draw(Graphics g, int[] origin, double scale) {;}
    public void initializeSpecies(Phase phase) {
        parentPhase = phase;
        firstAtom().r[0] = 0.5 - 0.5*diameter;
        lastAtom().r[0] = 0.5 + 0.5*diameter;
        firstAtom().r[1] = 0.0;
        lastAtom().r[1] = 0.0;
        firstAtom().setDiameter(1.0);
        lastAtom().setDiameter(1.0);
    }
  }

  class PotentialMomentumFlux extends PotentialHardDiskWall {
    
    public PotentialMomentumFlux() {
        super(0.0);
    }

    public double collisionTime(Atom atom1, Atom atom2) {
   
        Atom disk;
        AtomWall wall;
        if(atom2 instanceof AtomWall) {
           disk = atom1;
           wall = (AtomWall)atom2;
        }
        else {
           disk = atom2;
           wall = (AtomWall)atom1;
        }
        
        double time = Double.MAX_VALUE;
        int i = 0;
            double dr, t, dtdr;
            dr = wall.r[i] - disk.r[i];
            dtdr = 1.0/(disk.p[i]*disk.rm);
            t = dr*dtdr;

            if(t > 0.0) {time = Math.max(0.0,t-collisionDiameter/2.0*Math.abs(dtdr));}
            return time;
    }
    
    public void bump(Atom atom1, Atom atom2)  //this needs updating to check for isStationary
    {
        Atom disk;
        AtomWall wall;
        if(atom2 instanceof AtomWall) {
           disk = atom1;
//           wall = (AtomWall)atom2;
        }
        else {
           disk = atom2;
//           wall = (AtomWall)atom1;
        }
                
        momentumSum += Math.abs(disk.p[0]);
    }
        
  }        
}
