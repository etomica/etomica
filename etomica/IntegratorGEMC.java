package simulate;

import java.util.Random;

public class IntegratorGEMC extends Integrator {
    
    private transient final double[] dr = new double[Space.D];
    private final Random rand = new Random();
    public double maxRStep, maxVStep;
    public double pressure, betaMu, eBetaMu;
    private int freqDisplace, freqVolume, freqMolecule;
    private int iDisplace, iVolume, iMolecule, iTotal;
    public Phase secondPhase;
    
    public IntegratorGEMC() {
        super();
        freqDisplace = 1;
        freqVolume = freqMolecule = 1;
        maxRStep = 0.1;
        maxVStep = 0.1;
        nPhasesMax = 2;
        phase = new Phase[nPhasesMax];
    }
    
    public void doStep(double dummy) {
        int i = (int)(rand.nextDouble()*iTotal);
        if(i < iDisplace) {
            if(rand.nextDouble() < 0.5) {trialDisplace(firstPhase);}
            else {trialDisplace(secondPhase);}
        }
        else if(i < iVolume) {
            trialVolume();
        }
        else {
            if(rand.nextDouble() < 0.5) {trialExchange(firstPhase.firstSpecies,secondPhase.firstSpecies);}
            else {trialExchange(secondPhase.firstSpecies,firstPhase.firstSpecies);}
        }
    }
        
    private void trialDisplace(Phase p) {
        if(p.nMoleculeTotal == 0) {return;}
        double uOld, uNew;
        int i = (int)(rand.nextDouble()*p.nAtomTotal);
        Atom a = p.firstAtom();
        for(int j=i; --j>=0; ) {a = a.getNextAtom();}  //get ith atom in list
        uOld = a.potentialEnergy();
        Space.randomVector(dr, maxRStep, rand);
        a.displace(dr);
        uNew = a.potentialEnergy();
        if(uNew < uOld) {return;}   //accept
        if(uNew >= Double.MAX_VALUE || 
           Math.exp(-(uNew-uOld)/temperature) < rand.nextDouble()) { //reject
             a.replace();
             return;
        }           
    }
    
    private void trialVolume() {
        double v1Old = firstPhase.space.volume;
        double v2Old = secondPhase.space.volume;
        double hOld = firstPhase.potentialEnergy() + secondPhase.potentialEnergy();
        double vStep = (2.*rand.nextDouble()-1.)*maxVStep;
        double v1New = v1Old + vStep;
        double v2New = v2Old - vStep;
        double v1Scale = v1New/v1Old;
        double v2Scale = v2New/v2Old;
        double r1Scale = Math.pow(v1Scale,1.0/(double)Space.D);
        double r2Scale = Math.pow(v2Scale,1.0/(double)Space.D);
        firstPhase.space.inflate(r1Scale);
        secondPhase.space.inflate(r2Scale);
        for(Molecule m=firstPhase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
            m.inflate(r1Scale);
        }
        for(Molecule m=secondPhase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
            m.inflate(r2Scale);
        }
        double hNew = firstPhase.potentialEnergy() + secondPhase.potentialEnergy();
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/temperature+
                       firstPhase.nMoleculeTotal*Math.log(v1Scale) +
                       secondPhase.nMoleculeTotal*Math.log(v2Scale))
                < rand.nextDouble()) 
            {  //reject
              firstPhase.space.inflate(1.0/r1Scale);
              for(Molecule m=firstPhase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
                m.replace();
              }
              secondPhase.space.inflate(1.0/r2Scale);
              for(Molecule m=secondPhase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
                m.replace();
              }
            }
    }
    
    private void trialExchange(Species iSpecies, Species dSpecies) {
        double uNew, uOld;
        
        if(dSpecies.nMolecules == 0) {return;}
 
        iSpecies.parentPhase.space.randomVector(dr, rand);  //random point in volume
        Molecule iM = new Molecule(iSpecies,iSpecies.nAtomsPerMolecule);
        iM.translate(dr);
        uNew = iM.potentialEnergy();
        if(uNew == Double.MAX_VALUE) {return;}        //overlap
        
        int i = (int)(rand.nextDouble()*dSpecies.nMolecules);
        Molecule dM = dSpecies.firstMolecule;
        for(int j=i; --j>=0; ) {dM = dM.getNextMolecule();}
        uOld = dM.potentialEnergy();
        
        double bFactor = dSpecies.nMolecules/dSpecies.parentPhase.space.volume
                         * iSpecies.parentPhase.space.volume/(iSpecies.nMolecules+1)
                         * Math.exp(-(uNew-uOld)/temperature);
        if(bFactor > 1.0 || bFactor > rand.nextDouble()) {  //accept
            dSpecies.deleteMolecule(dM);
            iSpecies.addMolecule(iM);
        }                   
    }
            
    public final int getFreqVolume() {return freqVolume;}
    public final void setFreqVolume(int f) {freqVolume = f;}
    public final int getFreqMolecule() {return freqMolecule;}
    public final void setFreqMolecule(int f) {freqMolecule = f;}
    public final int getFreqDisplace() {return freqDisplace;}
    public final void setFreqDisplace(int f) {freqDisplace = f;}
        
    public final double getMaxRStep() {return maxRStep;}
    public final void setMaxRStep(double s) {maxRStep = s;}
    
    public final double getMaxVStep() {return maxVStep;}
    public final void setMaxVStep(double s) {maxVStep = s;}
    
    public void initialize() {
        iDisplace = freqDisplace * firstPhase.nMoleculeTotal;
        iVolume = iDisplace + freqVolume;
        iMolecule = iVolume + freqMolecule*firstPhase.nMoleculeTotal;    
        iTotal = iMolecule;
    }

}