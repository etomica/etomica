package etomica.data.meter;

import etomica.AtomPair;
import etomica.Default;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.MeterAbstract;
import etomica.Phase;
import etomica.Space;
import etomica.atom.iterator.ApiLeafAtoms;
import etomica.atom.iterator.AtomsetIteratorPhaseDependent;
import etomica.math.SphericalHarmonics;
import etomica.space.CoordinatePair;
import etomica.space.Vector;
import etomica.units.Dimension;

 /** The Bond Order Parameter Ql provides a metric that indicates the crystallinity of a phase.
   * Appropriate for 3-dimensional system only.
   * Refer: Journal of Chemical Physics Vol.104 No.24,22nd June,1996__ Rate of crystal nucleation
   *
   * @author Jhumpa Adhikari
   */

public class MeterBondOrderParameterQ extends MeterAbstract implements EtomicaElement {
	
    public MeterBondOrderParameterQ(Space space) {
        super(1);
        setL(6);
        setR2Cut(Math.pow(5.0*Default.ATOM_SIZE, 2));
        setLabel("Bond Q Order Parameter");
        cPair = space.makeCoordinatePair();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Order parameter used to identify crystallinity of phases (3-D only)");
        return info;
    }

    /**
     * Returns the value of the bond-order parameter for the given phase
     * in its current configuration.  Returned array has only one element.
     */
    public double[] getData(Phase phase) {
        int nbSum = 0;
        for(int m=-L; m<=L; m++) {
            int idx = m+L;
            Qreal[idx] = 0.0;
            Qimag[idx] = 0.0;
        }
        cPair.setNearestImageTransformer(phase.boundary());
        pairIterator.setPhase(phase);
        pairIterator.reset();
        while(pairIterator.hasNext()) {
            AtomPair pair = (AtomPair)pairIterator.next();
        	cPair.reset(pair.atom0.coord,pair.atom1.coord);
        	double r2 = cPair.r2();
            if(r2 < r2Cut) {
                nbSum += 2;
                Vector rVec = cPair.dr();
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
        value[0] = Math.sqrt(coeff*QL)/(double)nbSum;
        return value;
        
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
        
    public Dimension getDimension() {return Dimension.NULL;}

    /**
     * Sets the iterator that gives the atoms over which the order parameter
     * will be calculated.  Default iterator is ApiLeafAtoms.
     * @param iter
     */
    public void setIterator(AtomsetIteratorPhaseDependent iter) {
    	if(iter.nBody() != 2) throw new IllegalArgumentException("Illegal attempt to use a non-pair iterator");
    	pairIterator = iter;
    }
    /**
     * @return the iterator giving the atoms over which the order parameter
     * is defined.
     */
    public AtomsetIteratorPhaseDependent getIterator() {
    	return pairIterator;
    }
    
    public double getR2Cut(){
    	return r2Cut;
    }
    
    public void setR2Cut(double r2c){
    	r2Cut = r2c;
    }
	
    private SphericalHarmonics sh;
    private double[] Qreal, Qimag;
    private int L;
    private AtomsetIteratorPhaseDependent pairIterator = new ApiLeafAtoms();
    private double r2Cut;
    private double[] rThetaPhi = new double[3];
    private final double[] value = new double[1];
    private double coeff;
    private final CoordinatePair cPair;
    

/*    public static void main(String[] args) {
        
        Default.ATOM_SIZE = 1.0;
        Simulation sim = new Simulation(new Space3D());
        Simulation.instance = sim;
        
        Species species = new SpeciesSpheres(32);
        P2HardSphere potential = new P2HardSphere();
        Integrator integrator = new IntegratorHard();
        Controller controller = new Controller();
        Meter meter = new MeterBondOrderParameterQ();
        DisplayBox box = new DisplayBox();
        Phase phase = new Phase();
        Configuration configuration = new ConfigurationFcc(sim.space());
        DisplayPhase display = new DisplayPhase();
        
        phase.setConfiguration(configuration);
        
        box.setMeter(meter);
        box.setPrecision(8);
        
        sim.elementCoordinator.go();
        phase.setDensity(1.1);
        meter.updateSums();
        box.doUpdate();
        display.setScale(2.0);

        Simulation.makeAndDisplayFrame(sim);
    }//end of main
    */
}//end of MeterBondOrderParameterQ

