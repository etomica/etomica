package etomica;

public class MCMoveBiasUB extends MCMove {
    
    private Space.Vector e1 =  Simulation.instance.space().makeVector();
    private  Space.Vector e2 =  Simulation.instance.space().makeVector();
    private double wellCutoffSquared;
    private BiasVolumeUB biasVolume;
    private int Ni, Nai;
    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
    private final PotentialCalculation.EnergySum energy = new PotentialCalculation.EnergySum();

    public MCMoveBiasUB(BiasVolumeUB bv) {
        super();
        setStepSizeMax(Default.BOX_SIZE);
        setStepSizeMin(0.0);
        setStepSize(Default.ATOM_SIZE);
        setPerParticleFrequency(true);
        biasVolume = bv;
    }
    
   public void setPhase(Phase p) {
        phase = p;
        bv.setPhase(p);
    }
    
    public void thisTrial() {
        
        if(phase.moleculeCount() < 2) return;
        /*
        choose bonding or unbonding
        */
        double d = rand.nextDouble();
        double r2;
     
        double uOld, uNew;
         
        if (d < 0.5) { // bonding
         
            Atom atom = phase.randomMolecule();
            uOld = potential.set(phase).calculate(iteratorDirective.set(atom), energy.reset()).sum();
            Atom a1 =bv.biasInsert();
         
        if (a1 != null ) {
            
            Ni =  bv.getNi();
            Nai= bv.getNai();
            if ( acceptanceParameter(a1, Ni, Nai, uOld) < rand.nextDouble()) // reject
            {
              
             
                a1.replace();
                phase.iteratorFactory().moveNotify(a1);
                return;
            }
            nAccept++;   //accept 
        }
        }//end unbonding
        else { // bonding
         
        Atom a1 = bv.biasRemove();
        Ni =  bv.getNi();
        Nai= bv.getNai();
            
             
            if ( a1 != null) {  
            
            
            
            if ( acceptanceParameter(a1, Ni, Nai, uOld) < rand.nextDouble()) //reject
            {
              
                a1.replace();
                phase.iteratorFactory().moveNotify(a1);
                return;
            }
            
            nAccept++;   //accept
            }
        }
     
     
      
    }
    
 
      /* Method for acceptance parameter */
      double acceptanceParameter(Atom a, double Ni, double Nai, double uOld){
        
       
        double N = phase.atomCount();
        double Naj, Nj,phi; 
        double uNew = phase.energy.meterPotential().currentValue(a);
        phi = bv.volumeFraction();
        
        Naj=bv.numberOfBondedMolecules();
        Nj =bv.returnBondingRegionMoleculeCount(a);
     
        
       
       return (((N-1)*phi*(rand.nextDouble())/Naj + Ni)/((N-1.0)*phi*(rand.nextDouble())/Nai + Nj))*
             Math.exp(-(uNew-uOld)/parentIntegrator.temperature);    
       
       
       
        
        
        
      }
}