package etomica;

import etomica.units.Dimension;
import etomica.*;
import Jama.util.*;
import Jama.*;

public class MCMoveDeform extends MCMove {
    
    protected double pressure;
    public Space.Tensor pressureTensor; 
    public Space.Tensor dualTensor;
    public Space.Tensor stressTensor;
    
    public Space.Tensor deformationTensor;
    public Space.Tensor eOld;
    public Space.Tensor eNew;
    private Space.Tensor strain;

    public Space.Tensor HOld;
    public Space.Tensor HNew;
    public boolean firstTime=true;
    public int accept=1;
    public int reject=1;
    private double hOld;
    
    public double omega0;
    protected PhaseActionSlantStress phaseActionSlantStress;
    private final IteratorDirective iteratorDirective = new IteratorDirective();
    private AtomIterator affectedAtomIterator;

    public MCMoveDeform(IntegratorMC parentIntegrator) {
        super(parentIntegrator);
        pressureTensor = parentIntegrator.parentSimulation().space().makeTensor();
        dualTensor = parentIntegrator.parentSimulation().space().makeTensor();
            stressTensor = parentIntegrator.parentSimulation().space().makeTensor();
            deformationTensor = parentIntegrator.parentSimulation().space().makeTensor();
        eOld = parentIntegrator.parentSimulation().space().makeTensor();
        eNew = parentIntegrator.parentSimulation().space().makeTensor();
        strain = parentIntegrator.parentSimulation().space().makeTensor();
        HOld = parentIntegrator.parentSimulation().space().makeTensor();
        HNew = parentIntegrator.parentSimulation().space().makeTensor();

        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.1);
        setPressure(Default.PRESSURE);
        setPressureTensor(Default.PRESSURE);
        setStressTensor(dualTensor);
            HOld.E(0.0);
            HNew.E(0.0);
                for(int i=0;i<=2;i++){
                HOld.setComponent(i, i, Default.BOX_SIZE);
                HNew.setComponent(i, i, Default.BOX_SIZE);
                } 
            defineStrains(HOld, HNew);
   }
    
    public void setPhase(Phase p) {
        if(p == null) return;
        super.setPhase(p);
        phaseActionSlantStress = new PhaseActionSlantStress(phase);
        affectedAtomIterator = phase.makeAtomIterator();
        setInitialVolume(p);
    }
        
    public double trace(Space.Tensor st, Space.Tensor p, Space.Tensor ep){
      st.ME(p); 
//         for(int i=0;i<3;i++){
//         for(int j=0;j<3;j++){
//            System.out.println(st.component(i,j)+" "+p.component(i,j));
//         }         
//         }
         
    return  st.component(0,0)*ep.component(0,0)+st.component(0,1)*ep.component(1,0)+st.component(0,2)*ep.component(2,0)+
            st.component(1,0)*ep.component(0,1)+st.component(1,1)*ep.component(1,1)+st.component(1,2)*ep.component(2,1)+
            st.component(2,0)*ep.component(0,2)+st.component(2,1)*ep.component(1,2)+st.component(2,2)*ep.component(2,2);
    }
    
    public void defineStrains (Space.Tensor ho, Space.Tensor h){
        
        double[][] vals0 = {{ho.component(0, 0),ho.component(0, 1),ho.component(0, 2)}
                           ,{ho.component(1, 0),ho.component(1, 1),ho.component(1, 2)}
                           ,{ho.component(2, 0),ho.component(2, 1),ho.component(2, 2)}};
        double[][] vals1 = {{h.component(0, 0),h.component(0, 1),h.component(0, 2)}
                           ,{h.component(1, 0),h.component(1, 1),h.component(1, 2)}
                           ,{h.component(2, 0),h.component(2, 1),h.component(2, 2)}};
      
        Matrix Ho = new Matrix(vals0);
        Matrix H = new Matrix(vals1);
        Matrix G = H.transpose().times(H);
        Matrix Hoit=Ho.transpose().inverse();
        Matrix Hoi=Ho.inverse();
        Matrix HoitGHoi = Hoit.times(G).times(Hoi);
      if(firstTime){
            eOld.E(0.0);                            
            firstTime=false;        
      }else{
            eNew.setComponent(0,0, 0.5*(HoitGHoi.get(0,0)-1.0));
            eNew.setComponent(0,1, 0.5*(HoitGHoi.get(0,1)));
            eNew.setComponent(0,2, 0.5*(HoitGHoi.get(0,2)));
            eNew.setComponent(1,0, 0.5*(HoitGHoi.get(1,0)));
            eNew.setComponent(1,1, 0.5*(HoitGHoi.get(1,1)-1.0));
            eNew.setComponent(1,2, 0.5*(HoitGHoi.get(1,2)));
            eNew.setComponent(2,0, 0.5*(HoitGHoi.get(2,0)));
            eNew.setComponent(2,1, 0.5*(HoitGHoi.get(2,1)));
            eNew.setComponent(2,2, 0.5*(HoitGHoi.get(2,2)-1.0)); 
    }
    }

    public boolean doTrial() {
        setStressTensor(dualTensor);
        double omega=phase.volume();
  
        double uOld = potential.set(phase).calculate(iteratorDirective, energy.reset()).sum();
        hOld = uOld +omega0*trace(stressTensor, pressureTensor, eOld);//+ pressure*(omega-omega0)
        setStressTensor(dualTensor);

        double rScaleN1 = (2.*Math.random()-1.)*stepSize+1;  
        
        phaseActionSlantStress.setScale(1 , rScaleN1);
//        phaseActionSlantStress.setScale(3 , rScaleN1);
        phaseActionSlantStress.attempt();
        return true;
    }
    
    public double lnProbabilityRatio() {
        HNew.setComponent(1,0, (((Space.Boundary.Deformable)phase.boundary()).getDeformationTensor().component(0,1)));  
//        System.out.println(" Component of H (1,0)"+HNew.component(1,0));
//        HNew.setComponent(0,1, (((Space.Boundary.Deformable)phase.boundary()).getDeformationTensor().component(1,0)));  
        defineStrains(HOld, HNew);          

        double uNew = potential.set(phase).calculate(iteratorDirective, energy.reset()).sum();
        double hNew = uNew +omega0*trace(stressTensor, pressureTensor, eNew);//+ pressure*(omega-omega0)
      
//  if(hNew!=0.0){System.out.println(" Energy "+hOld+" "+hNew);}
     
        return -(hNew-hOld)/parentIntegrator.temperature;
    }
    
    public double lnTrialRatio() {return 0.0;}
    
    public void rejectNotify() {
              phaseActionSlantStress.undo();
    }
    
    public void acceptNotify() {
            eOld.E(eNew);
            HOld.E(HNew);
            strain.E(eNew);   
    }

    public Space.Tensor getStrain () {
       return strain;
    }
    
    public AtomIterator affectedAtoms(Phase phase) {
        throw new RuntimeException("affectedAtoms not implemented in MCMoveDeform");
    }

    public void setPressure(double p) {pressure = p;}
    public void setPressureTensor(double p) {pressureTensor.setComponent(0,0, p);
                                             pressureTensor.setComponent(1,1, p);
                                             pressureTensor.setComponent(2,2, p);
    }
    public final double getPressure() {return pressure;}
    public final Space.Tensor getPressureTensor() {return pressureTensor;}
    public final double pressure() {return pressure;}
    public Dimension getPressureDimension() {return Dimension.PRESSURE;}
    public final void setLogPressure(int lp) {setPressure(Math.pow(10.,(double)lp));}
    
    public Space.Tensor dualTensor() {return dualTensor;}

    public void setStressTensor(Space.Tensor s) { 
         stressTensor.setComponent(0,0,(s.component(0,0)));//Pxx
         stressTensor.setComponent(0,1,(s.component(0,1)));
         stressTensor.setComponent(0,2,(s.component(0,2)));
         stressTensor.setComponent(1,0,(s.component(1,0)));
         stressTensor.setComponent(1,1,(s.component(1,1)));//Pyy
         stressTensor.setComponent(1,2,(s.component(1,2)));
         stressTensor.setComponent(2,0,(s.component(2,0)));
         stressTensor.setComponent(2,1,(s.component(2,1)));
         stressTensor.setComponent(2,2,(s.component(2,2)));//Pzz
        }
    public Space.Tensor stressTensor() {return stressTensor;}
    
    public void setDeformationTensor(Space.Tensor s) { deformationTensor.E(s);}    
    public Space.Tensor deformationTensor() {return deformationTensor;}
    
    public void setInitialVolume(Phase p) { omega0=p.volume(); 
    System.out.println("volume "+omega0);
    }

}