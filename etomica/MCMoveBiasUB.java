package etomica;
import etomica.association.*;


public class MCMoveBiasUB extends etomica.MCMove {
    
    private Space.Vector e1 =  Simulation.instance.space().makeVector();
    private  Space.Vector e2 =  Simulation.instance.space().makeVector();
    private etomica.association.BiasVolume bv;
    private AssociationManager associationManager;
    private int Ni, Nai;
    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private final double thetac = etomica.units.Degree.UNIT.toSim(27.0);
    private MCMoveEvent evt;
    
    private Atom atomA;
   
    public MCMoveBiasUB(IntegratorMC parentIntegrator, etomica.association.BiasVolume biasVolume, Phase phase) {
        super(parentIntegrator);
        setStepSizeMax(Default.BOX_SIZE);
        setStepSizeMin(0.0);
        setStepSize(Default.ATOM_SIZE);
        setPerParticleFrequency(true);
        setPhase(phase);
        atomA = (Atom)phase.firstAtom();
       
        bv=biasVolume;
       // associationManager = new AssociationManager(affectedAtomIterator, bv);
       //evt = new MCMoveEvent(this);
        
    }
    public void setAssociationManager(etomica.association.AssociationManager am){
        associationManager= am;   
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
         
        if (d < 0.5) { // bonding
          System.out.println(" I am in bonding");
            atomA = phase.randomMolecule();
            Atom atomB = atomA;
            while(atomB == atomA) atomB = phase.randomMolecule();
            
            uOld = potential.set(phase).calculate(iteratorDirective.set(atomA), energy.reset()).sum();

            ni = associationManager.associationCount(atomA);
            Nai = associationManager.associatedAtomCount();
            deltai = (ni == 0) ? 0 : 1;

            bv.biasInsert(atomA, atomB);
            
            uNew = potential.set(phase).calculate(iteratorDirective.set(atomA), energy.reset()).sum();
            if(uNew == Double.MAX_VALUE) {
                //evt.acceptedMove=false;        
                System.out.println(" yo");
                atomA.coord.replace(); 
                return false;
            }
            
            deltaU = uNew - uOld;
            System.out.println(" deltaU" + deltaU);
            //evt.acceptedMove=true;
            //associationManager.mcMoveAction(evt);
            //associationManager.update(atomA);
            nj = associationManager.associationCount(atomA);
            Naj = associationManager.associatedAtomCount();
            deltaj = (nj == 0) ? 0 : 1;
            //return true;
        }//end bonding
        else { // unbonding
       
            atomA = associationManager.randomAssociatedAtom();
            
            if(atomA != null) {
                 System.out.println(" I am in unbonding");
            ni = associationManager.associationCount(atomA);
            Nai = associationManager.associatedAtomCount();
            deltai = (ni == 0) ? 0 : 1;
            uOld = potential.set(phase).calculate(iteratorDirective.set(atomA), energy.reset()).sum();

            atomA.coord.translateTo(phase.randomPosition());
           // System.out.println(" atomA x" +((Space2D.Vector)atomA.coord.position()).component(0));
            // System.out.println(" atomA y" +((Space2D.Vector)atomA.coord.position()).component(1));
              e1.E(0.);
             e1.setComponent(0,1);
            ((Space.Coordinate.Angular)atomA.coord).orientation().convertToSpaceFrame(e1);
      
            e1.randomRotate(thetac);
         
            
            uNew = potential.set(phase).calculate(iteratorDirective.set(atomA), energy.reset()).sum();
            if(uNew == Double.MAX_VALUE) {
                
                atomA.coord.replace(); 
                return false;
            }
            
            deltaU = uNew - uOld;
             //evt.acceptedMove=true;
            //associationManager.mcMoveAction(evt);
            nj = associationManager.associationCount(atomA);
            Naj = associationManager.associatedAtomCount();
            deltaj = (nj == 0) ? 0 : 1;
            }
            else return false;
        }//end unbonding
      
        double phi = bv.biasVolume()/phase.volume();
        double chi;
        
        if(Naj == 0) chi = ni*Nai/((N-1)*phi) * Math.exp(-deltaU/parentIntegrator.temperature);
        else if(Nai == 0) chi = (N-1)*phi/(nj*Naj) * Math.exp(-deltaU/parentIntegrator.temperature);
        else chi = ((N-1)*phi*deltaj/Naj + ni)/((N-1)*phi*deltai/Nai + nj) *
                      Math.exp(-deltaU/parentIntegrator.temperature);    

        if ( chi < Simulation.random.nextDouble()) {// reject
          System.out.println( "inside chi");
            atomA.coord.replace();
             //evt.acceptedMove=false;
            //associationManager.mcMoveAction(evt);
            //associationManager.update(atomA);
            return false;
        }
        System.out.println( "outside chi");
        
        //nAccept++;   //accept 
        return true;
    }//end thisTrial
    
  public final AtomIterator affectedAtoms() {
        affectedAtomIterator.setBasis(atomA);
        return affectedAtomIterator;
    }

}//end of MCMoveBiasUB