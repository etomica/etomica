package simulate;
import java.beans.*;

public class IntegratorVelocityVerlet extends Integrator {

    private transient final double dr[] = new double[Space.D];
    private transient final double dp[] = new double[Space.D];

    public IntegratorVelocityVerlet() {
        super();
        setNeighborListUpdateInterval(1); 
    }

    public void doStep(double tStep) {    
        double halfTStep = 0.5*tStep;
        for(int j=phase.speciesVector.size(); --j>=0; ) {
            Species species = (Species)phase.speciesVector.elementAt(j);
            if(species instanceof SpeciesDisk) {
                SpeciesDisk speciesDisk = (SpeciesDisk)species;
                double a1 = speciesDisk.rm*tStep;    // rm*tStep
                double a2 = a1*halfTStep;            // 0.5*rm*tStep^2
                for (int i=speciesDisk.nElements; --i>=0; ) {
                   MoleculeAtomic m = speciesDisk.molecule[i];
                   Space.uEa1Tv1Pa2Tv2(dr, a1, m.p, a2, m.f);
                   Space.uEa1Tv1(dp, halfTStep, m.f);
                   m.translate(dr);
                   m.accelerate(dp);
                }
            }
        }
        phase.updatedForces = false;
        phase.updatedKineticEnergy = false;
        phase.updateForces(); 
        for(int j = 0; j < phase.speciesVector.size(); j++) {
            Species species = (Species)phase.speciesVector.elementAt(j);
            if(species instanceof SpeciesDisk) {
                SpeciesDisk speciesDisk = (SpeciesDisk)species;
                double rm = speciesDisk.rm;
                for (int i = 0; i < speciesDisk.nElements; i++) {
                   MoleculeAtomic m = speciesDisk.molecule[i];
                   Space.uEa1Tv1(dp, halfTStep, m.f);
                   m.accelerate(dp);
                }
            }
        }
        phase.updatedKineticEnergy = false;
    } //end of doStep

    public void initialize() {
        phase.updateNeighbors();
        phase.updateForces();
    }

} //end of class