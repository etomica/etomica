package etomica;

import etomica.units.*;
import etomica.utility.*;

 /** The Bond Order Parameter Ql provides a metric that indicates the crystallinity of a phase.
   * Appropriate for 3-dimensional system only.
   * Refer:journal Of Chemical Physics Vol.104 No.24,22nd June,1996__ Rate of crystal nucleation
   *
   * @author Jhumpa Adhikari
   */

public class MeterBondOrderParameterQ extends Meter implements EtomicaElement 
{
    private SphericalHarmonics sh;
    private double[] Qreal, Qimag;
    private int L;
    private AtomPairIterator pairIterator;
    private double r2Cut;
    private double[] rThetaPhi = new double[3];
    private double coeff;
    
    public MeterBondOrderParameterQ() {this(Simulation.instance);}
    
    public MeterBondOrderParameterQ(Simulation sim) {
        super(sim);
        setL(6);
        setR2Cut(Math.pow(5.0*Default.ATOM_SIZE, 2));
        setLabel("Bond Q Order Parameter");
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Order parameter used to identify crystallinity of phases (3-D only)");
        return info;
    }

    public double currentValue() {
        int nbSum = 0;
        for(int m=-L; m<=L; m++) {
            int idx = m+L;
            Qreal[idx] = 0.0;
            Qimag[idx] = 0.0;
        }
        pairIterator.reset();
        while(pairIterator.hasNext()) {
            AtomPair pair = pairIterator.next();
            double r2 = pair.r2();
            if(r2 < r2Cut) {
                nbSum += 2;
                Space.Vector rVec = pair.dr();
                rVec.sphericalCoordinates(rThetaPhi);
                double theta = rThetaPhi[1];
                double phi = rThetaPhi[2];
                for(int m=-L; m<=L; m++) {
                    int idx = m+L;
                    double thetaC = Math.PI - theta;
                    double phiC = Math.PI + phi;
                    Qreal[idx] += SphericalHarmonics.realYm(L, m, theta, phi);
                    Qimag[idx] += SphericalHarmonics.imaginaryYm(L, m, theta, phi);
                    Qreal[idx] += SphericalHarmonics.realYm(L, m, thetaC, phiC);
                    Qimag[idx] += SphericalHarmonics.imaginaryYm(L, m, thetaC, phiC);
                }//end for
            }//end if
        }//end while
        double QL = 0;
        for(int m=-L; m<=L; m++) {
            int idx = m+L;
            QL += Qreal[idx]*Qreal[idx] - Qimag[idx]*Qimag[idx];
        }
        QL = Math.sqrt(coeff*QL)/(double)nbSum;
        return QL;
        
    }//end of currentValue
    
    public int getL(){return L;}
    /**
     * Sets the value of the parameter l, to indicate if Q4 or Q6 is to be computed.
     * Input of any value other than 4 causes L to be set to 6.
     */
    public void setL(int L){
        if(L != 4) L = 6;
        this.L = L;
        Qreal = new double[2*L + 1];
        Qimag = new double[2*L + 1];
        coeff = 4*Math.PI/(2*L + 1);
    }
        
    public Phase getPhase() {return super.getPhase();}
    public void setPhase(Phase p) {
        super.setPhase(p);
        if(pairIterator == null) pairIterator = new AtomPairIterator(p);
    }
    public Dimension getDimension() {return Dimension.NULL;}

    public void setIterator(AtomPairIterator iter) {pairIterator = iter;}
    public AtomPairIterator getIterator() {return pairIterator;}
    
    public double getR2Cut(){return r2Cut;}
    public void setR2Cut(double r2c){r2Cut = r2c;}
	
    public static void main(String[] args) {
        
        Default.ATOM_SIZE = 1.0;
        Simulation sim = new Simulation(new Space3D());
        Simulation.instance = sim;
        
        Species species = new SpeciesDisks(32);
        P2HardSphere potential = new P2HardSphere();
        Integrator integrator = new IntegratorHard();
        Controller controller = new Controller();
        Meter meter = new MeterBondOrderParameterQ();
        DisplayBox box = new DisplayBox();
        Phase phase = new Phase();
        Configuration configuration = new ConfigurationFcc();
        DisplayPhase display = new DisplayPhase();
        
        phase.setConfiguration(configuration);
        
        box.setMeter(meter);
        box.setPrecision(8);
        
        sim.elementCoordinator.go();
        phase.setDensity(1.1);
        meter.updateSums();
        box.doUpdate();
        display.setScale(2.0);
        Potential2.Agent potentialAgent = (Potential2.Agent)potential.getAgent(phase);
        potentialAgent.setIterator(new AtomPairIterator(phase));

        Simulation.makeAndDisplayFrame(sim);
    }//end of main
}

