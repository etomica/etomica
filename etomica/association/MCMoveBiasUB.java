package etomica.association;
import etomica.*;

public class MCMoveBiasUB extends MCMove {
    
    private Space.Vector e1 =  Simulation.instance.space().makeVector();
    private  Space.Vector e2 =  Simulation.instance.space().makeVector();
    private BiasVolume biasVolume;
    private AssociationManager associationManager;
    private int Ni, Nai;
    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();

    public MCMoveBiasUB(IntegratorMC parentIntegrator, BiasVolume bv) {
        super(parentIntegrator);
        setStepSizeMax(Default.BOX_SIZE);
        setStepSizeMin(0.0);
        setStepSize(Default.ATOM_SIZE);
        setPerParticleFrequency(true);
        biasVolume = bv;
        associationManager = new AssociationManager(bv);
    }
    
    public boolean thisTrial() {
        
        int N = phase.moleculeCount();
        if(N < 2) return false;
        /*
        choose bonding or unbonding
        */
        double d = Simulation.random.nextDouble();
        double r2;
        
        int ni, nj, deltai, deltaj, Naj;
        double deltaU;
     
        double uOld, uNew;
        Atom atomA;
         
        if (d < 0.5) { // bonding
         
            atomA = phase.randomMolecule();
            Atom atomB = atomA;
            while(atomB == atomA) atomB = phase.randomMolecule();
            
            uOld = potential.set(phase).calculate(iteratorDirective.set(atomA), energy.reset()).sum();

            ni = associationManager.associationCount(atomA);
            Nai = associationManager.associatedAtomCount();
            deltai = (ni == 0) ? 0 : 1;

            biasVolume.biasInsert(atomA, atomB);
         
            uNew = potential.set(phase).calculate(iteratorDirective.set(atomA), energy.reset()).sum();
            if(uNew == Double.MAX_VALUE) {
                atomA.coord.replace(); 
                return false;
            }
            
            deltaU = uNew - uOld;
            
            associationManager.update(atomA);
            nj = associationManager.associationCount(atomA);
            Naj = associationManager.associatedAtomCount();
            deltaj = (nj == 0) ? 0 : 1;
        }//end bonding
        else { // unbonding
        
            atomA = associationManager.randomAssociatedAtom();
            ni = associationManager.associationCount(atomA);
            Nai = associationManager.associatedAtomCount();
            deltai = (ni == 0) ? 0 : 1;
            uOld = potential.set(phase).calculate(iteratorDirective.set(atomA), energy.reset()).sum();

            atomA.coord.translateTo(phase.randomPosition());
            atomA.coord.rotateToRandom();
            
            uNew = potential.set(phase).calculate(iteratorDirective.set(atomA), energy.reset()).sum();
            if(uNew == Double.MAX_VALUE) {
                atomA.coord.replace(); 
                return false;
            }
            
            deltaU = uNew - uOld;
            
            associationManager.update(atomA);
            nj = associationManager.associationCount(atomA);
            Naj = associationManager.associatedAtomCount();
            deltaj = (nj == 0) ? 0 : 1;
        }//end unbonding
     
        double phi = biasVolume.biasVolume()/phase.volume();
        double chi;
        
        if(Naj == 0) chi = ni*Nai/((N-1)*phi) * Math.exp(-deltaU/parentIntegrator.temperature());
        else if(Nai == 0) chi = (N-1)*phi/(nj*Naj) * Math.exp(-deltaU/parentIntegrator.temperature());
        else chi = ((N-1)*phi*deltaj/Naj + ni)/((N-1)*phi*deltai/Nai + nj) *
                      Math.exp(-deltaU/parentIntegrator.temperature());    

        if ( chi < Simulation.random.nextDouble()) {// reject
            atomA.coord.replace();
            associationManager.update(atomA);
            return false;
        }
        nAccept++;   //accept
        return true;
    }//end thisTrial
 
}//end of MCMoveBiasUB