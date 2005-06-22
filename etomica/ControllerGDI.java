package etomica ;
import etomica.data.AccumulatorAverage;
import etomica.data.DataSourceScalar;
import etomica.modifier.ModifierGeneral;
import etomica.utility.Function;

/**
 * Controller that implements a Gibbs-Duhem integration for tracing a line of phase coexistence
 *
 * @author Jhumpa Adhikari
 */
 
public class ControllerGDI extends Controller implements EtomicaElement {
    
    private  double          independentVariable;   //Independent Variable
    private  double          dependentVariableOld ; //Dependent Variable at previous step
    private  double          dependentVariableNew ; //Dependent Variable at current step
    private  double          numerator;
    private  double          denominator;
    private  double          f;//slope of the coexistence line
    private  double          f0;//slope for Predictor step
    private  double          f1;//slope for Corrector step
    private  double          h = 0.0001;// step size
    private  double          factor;
    private  Integrator      integrator1,integrator2;
    private  ModifierGeneral       modifierIndependentVariable,modifierDependentVariable;
    private  DataSourceScalar     meter1phase1,meter1phase2,meter2phase1,meter2phase2;
    private  AccumulatorAverage accAve1phase1, accAve1phase2, accAve2phase1, accAve2phase2;
    private  DataManager AM1phase1, AM1phase2, AM2phase1, AM2phase2;
    private  double          i0 ;
    private  double          d0 ;
    private  double          independentVariableFinal ;//final value of independent variable
//    private  OutputFile      outputFile;
//    private  String          outputFileName;
    private int              noOfCorrectorIterations = 5;
    private int              pCycles = 1000;
    private int              cCycles = 1000;
    private int              prodCycles = 10000;
    private boolean          initialSlopeGiven = false;
    private Function iF = new Function.Identity();
    private Function dF = new Function.Identity();
    
    public ControllerGDI(){
        this(Simulation.instance);
    }
    public ControllerGDI(Simulation sim) {
        super(sim);
    }
    public ControllerGDI(Simulation sim, ModifierGeneral modifierIndependentVariable1, ModifierGeneral modifierDependentVariable1,
                            double i1, double d1, double finalvalue,
                            Integrator int1, Integrator int2,
                            DataSourceScalar m1, DataSourceScalar m2, DataSourceScalar m3, DataSourceScalar m4) {
        this(sim, modifierIndependentVariable1, modifierDependentVariable1,
                            i1,d1,finalvalue,
                            int1, int2,
                            m1, m2, m3, m4,
                            new Function.Identity(), new Function.Identity());
    }
    public ControllerGDI(Simulation sim, ModifierGeneral modifierIndependentVariable1, ModifierGeneral modifierDependentVariable1,
                            double i1, double d1, double finalvalue,
                            Integrator int1, Integrator int2,
                            DataSourceScalar m1, DataSourceScalar m2, DataSourceScalar m3, DataSourceScalar m4,
                            Function iFunction, Function dFunction) {
        super(sim);
        this.definition(modifierIndependentVariable1, modifierDependentVariable1,
                            i1,d1,finalvalue,
                            int1, int2,
                            m1, m2, m3, m4,
                            iFunction, dFunction);
    }
    
    public void definition(ModifierGeneral modifierIndependentVariable1, ModifierGeneral modifierDependentVariable1,
                            double i1,double d1,double finalvalue,
                            Integrator int1,Integrator int2,
                            DataSourceScalar m1, DataSourceScalar m2, DataSourceScalar m3, DataSourceScalar m4,
                            Function iFunction, Function dFunction) {
        iF = iFunction;
        dF = dFunction;
        modifierIndependentVariable = modifierIndependentVariable1;
        modifierDependentVariable   = modifierDependentVariable1;
        integrator1 = int1;
        integrator2 = int2;
        meter1phase1 = m1;
        meter1phase2 = m2;
        meter2phase1 = m3;
        meter2phase2 = m4;
        accAve1phase1 = new AccumulatorAverage();
        accAve1phase2 = new AccumulatorAverage();
        accAve2phase1 = new AccumulatorAverage();
        accAve2phase2 = new AccumulatorAverage();
        AM1phase1 = new DataManager(meter1phase1, new Accumulator[] {accAve1phase1});
        AM1phase2 = new DataManager(meter1phase2, new Accumulator[] {accAve1phase2});
        AM2phase1 = new DataManager(meter2phase1, new Accumulator[] {accAve2phase1});
        AM2phase2 = new DataManager(meter2phase2, new Accumulator[] {accAve2phase2});
        independentVariableFinal = iF.f(finalvalue);
        i0 = i1;
        d0 = d1;
        //should remove integrators, in case any were added previously, but don't have a remove method (yet)
        add(integrator1);
        add(integrator2);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Performs a Gibbs-Duhem integration phase-coexistence series");
        return info;
    }
    
    
    public void setCorrectorIterations(int steps){
        noOfCorrectorIterations = steps;
    }
    
    public int getCorrectorIterations(){
        return noOfCorrectorIterations;
    }
    
    public void setPredictorCycles(int c) {pCycles = c;}
    public int getPredictorCycles() {return pCycles;}
    public void setCorrectorCycles(int c) {cCycles = c;}
    public int getCorrectorCycles() {return cCycles;}
    public void setProductionCycles(int c) {prodCycles = c;}
    public int getProductionCycles() {return prodCycles;}
   
   public double getStepSize(){
        return h;
   }
   
   public void setStepSize(double steps){
        h = steps;
   }
   
   public void setInitialSlope(double f) {
        f0 = f;
        initialSlopeGiven = true;
   }
   public double getInitialSlope() {return f0;}
       
    private void doSimulation(int nCycles) {
        setMaxSteps(nCycles);
        meter1phase1.reset();
        meter1phase2.reset();
        meter2phase1.reset();
        meter2phase2.reset();
        currentEvent.setStarting(true);
        fireEvent(currentEvent);
        integrator1.start(); //start integrator on its own thread
        integrator1.join();  //wait until integrator finishes
        integrator2.start(); //don't start 2 until 1 is done, since thread safety is not assured in many classes
        integrator2.join();
        currentEvent.setStarting(false);
        fireEvent(currentEvent);
    }
    
    private double slope() {
        factor = dF.dfdx(modifierDependentVariable.getValue())/iF.dfdx(modifierIndependentVariable.getValue());
        numerator   = accAve1phase1.average()[0] - accAve1phase2.average()[0];
        denominator = accAve2phase1.average()[0] - accAve2phase2.average()[0];
        numerator *= factor;
        return numerator/denominator;
    }
    
    private void initial() {
//        outputFileName = "data4.txt";
//        outputFile = new OutputFile(outputFileName);
 
        independentVariable = iF.f(i0);
        modifierIndependentVariable.setValue(i0);
        modifierDependentVariable.setValue(d0);

        if(!initialSlopeGiven) {
            doSimulation(prodCycles);//simulate at initial conditions
            f0 = slope();
        }
        f1 = f0;
    }
    
    private void predictor() {
        f0 = f1;
        dependentVariableOld  = dF.f(modifierDependentVariable.getValue());
        dependentVariableNew = dependentVariableOld + f0*h;
        modifierDependentVariable.setValue(dF.inverse(dependentVariableNew));
        doSimulation(pCycles);
        f1 = slope();
    }
    
    
    private void corrector(int nCycles) {
        dependentVariableNew = dependentVariableOld + (f0 + f1)*0.5*h;
        modifierDependentVariable.setValue(dF.inverse(dependentVariableNew));
        doSimulation(nCycles);
        f1 = slope();
    }
    
//     private void output() {
//        outputFile.println(modifierIndependentVariable.getValue() +"\t"+ modifierDependentVariable.getValue());
//     }
 
     public void run() {
        
        initial();
        
        while(independentVariable < independentVariableFinal){
            
            independentVariable += h;
            modifierIndependentVariable.setValue(iF.inverse(independentVariable));
            
            currentEvent = predictorEvent;
            predictor();
            
            currentEvent = correctorEvent;
            for(int i = 1; i <= noOfCorrectorIterations; i++){
                corrector(cCycles);
            }
            currentEvent = productionEvent;
            corrector(prodCycles);
//            output();
        }
//        outputFile.close();
     }
     
     //Controller events used to signal the start/end of different simulation runs
     private Event currentEvent;
     private final Event predictorEvent = new Event(true, false, false);
     private final Event correctorEvent = new Event(false,true,false);
     private final Event productionEvent = new Event(false, false, true);
     /**
      * Event class that can indicate whether the controller is starting/finishing a predictor/corrector/production simulation.
      * Fired before and after performing a simulation for all phases.
      */
     public final class Event extends ControllerEvent {
        private final boolean predictor;
        private final boolean corrector;
        private final boolean production;
        public Event(boolean pred, boolean corr, boolean prod) {
            super(ControllerGDI.this);
            predictor = pred;
            corrector = corr;
            production = prod;
        }
        private boolean starting;
        public void setStarting(boolean b) {starting = b;}
        public boolean isCorrector() {return corrector;}
        public boolean isPredictor() {return predictor;}
        public boolean isProduction() {return production;}
        public boolean isStarting() {return starting;}
        public boolean isCompleted() {return !starting;}
     }
}