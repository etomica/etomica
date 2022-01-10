package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.*;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class SimLJPropsNormalModes extends Simulation {

    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Primitive primitive;
    protected MCMoveAtomCoupled atomMove;
    protected PotentialMasterList potentialMaster;
    protected final CoordinateDefinitionLeaf coordinateDefinition;
    public SpeciesGeneral species;
    protected  CoordinateDefinition.BasisCell[] cells0;


    public SimLJPropsNormalModes(Space _space, int numAtoms, double density, double temperature, double rC, boolean isLRC, double strain_x) {
        super(_space);
        setRandom(new RandomMersenneTwister(1)); // set seed
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        potentialMaster = new PotentialMasterList(this, space);

        double L = Math.pow(4.0/density, 1.0/3.0);
        System.out.println("rC = "+rC);
        int n = (int) Math.round(Math.pow(numAtoms/4, 1.0/3.0));

        nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        Vector[] cellDim = new Vector[3];
        cellDim[0] = Vector.of(new double[]{L, 0, 0});
        cellDim[1] = Vector.of(new double[]{0, L, 0});
        cellDim[2] = Vector.of(new double[]{0, 0, L});
        primitive = new PrimitiveGeneral(space, cellDim);
        boundary = new BoundaryDeformableLattice(primitive, nCells);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        CoordinateDefinitionLeaf coordinateDefinition0 = new CoordinateDefinitionLeaf(box, primitive, basisFCC, space);
        coordinateDefinition0.initializeCoordinates(new int[]{n, n, n});
        cells0 = coordinateDefinition0.getBasisCells();

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        atomMove = new MCMoveAtomCoupled(potentialMaster, meterPE, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(atomMove);


//        Potential2SoftSpherical potential = new P2LennardJones(space, 1.0, 1.0);
//        potential = new P2SoftSphericalTruncated(space, potential, rC);
        Potential2SoftSpherical potential = new P2SoftSphere(space);
        potential = new P2SoftSphericalTruncated(space, potential, rC);
        atomMove.setPotential(potential);

        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});

        potentialMaster.lrcMaster().setEnabled(isLRC);

        int cellRange = 4;
        potentialMaster.setRange(rC);
        potentialMaster.setCellRange(cellRange);
        potentialMaster.getNeighborManager(box).reset();
        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
        if (potentialCells < cellRange*2+1) {
            throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
        }
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addActivity(new ActivityIntegrate(integrator));

        if(!isLRC){
            ((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundary.getBoxSize().getX(0));
        }

        if (strain_x != 0) {
            boundary.setBoxSize(Vector.of(new double[]{(1+strain_x)*n*L, n*L, n*L}));
            cellDim[0] = Vector.of(new double[]{(1+strain_x)*L, 0, 0});
            cellDim[1] = Vector.of(new double[]{0, L, 0});
            cellDim[2] = Vector.of(new double[]{0, 0, L});
            primitive = new PrimitiveGeneral(space, cellDim);
            coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basisFCC, space);
            coordinateDefinition.initializeCoordinates(new int[]{n, n, n});

        } else {
            coordinateDefinition = coordinateDefinition0;
        }

    }

    public void initialize(long initSteps) {
        this.getController().runActivityBlocking(new ActivityIntegrate(this.integrator, initSteps));
        integrator.getMoveManager().setEquilibrating(false);
    }

    public static void main(String[] args) {
        SimOverlapParam params = new SimOverlapParam();
        ParseArgs.doParseArgs(params, args);

        double density = params.density;
        double rC = params.rC;
        long numSteps = params.numSteps;
        final int numAtoms = params.numAtoms;
        final int nC = params.nC;
        final int nBasis = params.nBasis;
        double temperature = params.temperature;
        double dbP = params.dbP;
        double dP = dbP*temperature;
        boolean isLRC = params.isLRC;
        double strain_x = params.strain_x;
        //Ulat=c2*rho^2 + c4*rho^4
        //at rC=3.0: 
        double c2=-14.195367616901942, c4=6.065858479622598;
        // Plat=-dU/dV = 
        // Plat = -dUlat/dV = + rho^2 dulat/drho = 2*c2*rho^3 + 4*c4*rho^5
        //dPlat/dV = -1/N* (6*c2*rho^4 + 20*c4*rho^6)
        //Blat = -V dPlat/dV
        //     = 1/rho* (6*c2*rho^4 + 20*c4*rho^6)
//        double dPLat =   1/density*( 6*c2*Math.pow(density,4.0)+20*c4*Math.pow(density, 6.0));

        String degK_k_g_eVecs = "degK_k_g_eVecs";
        File fileDegK_k_g_eVecs = new File(degK_k_g_eVecs);

        if(!fileDegK_k_g_eVecs.exists()){
            throw new RuntimeException("A file Does Not Exist!");
        }
        double[][] A;
        try{
            FileReader readerDegK_k_g_eVecs = new FileReader(degK_k_g_eVecs);
            BufferedReader buffDegK_k_g_eVecs = new BufferedReader(readerDegK_k_g_eVecs);
            String line;
            int nRows = 0;
            while((line = buffDegK_k_g_eVecs.readLine()) != null){
                nRows += 1;
            }
            int nCells = nC*nC*nC;
            int nColumns = 4+3*nBasis*(1+2*nBasis*3);
            A = new double[nRows][nColumns];
            int i = 0;
            readerDegK_k_g_eVecs = new FileReader(degK_k_g_eVecs);
            buffDegK_k_g_eVecs = new BufferedReader(readerDegK_k_g_eVecs);
            while((line = buffDegK_k_g_eVecs.readLine()) != null){
                line = line.trim();
                String[] v = line.split(" +");
                for (int j=0;j<nColumns;j++){
                    A[i][j] = Double.parseDouble(v[j]);
                }
                i+=1;
            }//while
        }catch(IOException ex){
            throw new RuntimeException(ex);
        }


        degK_k_g_eVecs = "degK_k_g_eVecs11";
        fileDegK_k_g_eVecs = new File(degK_k_g_eVecs);
        if(!fileDegK_k_g_eVecs.exists()){
            throw new RuntimeException("A file Does Not Exist!");
        }
        double[][] B;
        try{
            FileReader readerDegK_k_g_eVecs = new FileReader(degK_k_g_eVecs);
            BufferedReader buffDegK_k_g_eVecs = new BufferedReader(readerDegK_k_g_eVecs);
            String line;
            int nRows = 0;
            while((line = buffDegK_k_g_eVecs.readLine()) != null){
                nRows += 1;
            }
            int nCells = nC*nC*nC;
            int nColumns = 4+3*nBasis*(1+2*nBasis*3);
            B = new double[nRows][nColumns];
            int i = 0;
            readerDegK_k_g_eVecs = new FileReader(degK_k_g_eVecs);
            buffDegK_k_g_eVecs = new BufferedReader(readerDegK_k_g_eVecs);
            while((line = buffDegK_k_g_eVecs.readLine()) != null){
                line = line.trim();
                String[] v = line.split(" +");
                for (int j=0;j<nColumns;j++){
                    B[i][j] = Double.parseDouble(v[j]);
                }
                i+=1;
            }//while
        }catch(IOException ex){
            throw new RuntimeException(ex);
        }


        degK_k_g_eVecs = "degK_k_g_eVecs22";
        fileDegK_k_g_eVecs = new File(degK_k_g_eVecs);
        if(!fileDegK_k_g_eVecs.exists()){
            throw new RuntimeException("A file Does Not Exist!");
        }
        double[][] C;
        try{
            FileReader readerDegK_k_g_eVecs = new FileReader(degK_k_g_eVecs);
            BufferedReader buffDegK_k_g_eVecs = new BufferedReader(readerDegK_k_g_eVecs);
            String line;
            int nRows = 0;
            while((line = buffDegK_k_g_eVecs.readLine()) != null){
                nRows += 1;
            }
            int nCells = nC*nC*nC;
            int nColumns = 4+3*nBasis*(1+2*nBasis*3);
            C = new double[nRows][nColumns];
            int i = 0;
            readerDegK_k_g_eVecs = new FileReader(degK_k_g_eVecs);
            buffDegK_k_g_eVecs = new BufferedReader(readerDegK_k_g_eVecs);
            while((line = buffDegK_k_g_eVecs.readLine()) != null){
                line = line.trim();
                String[] v = line.split(" +");
                for (int j=0;j<nColumns;j++){
                    C[i][j] = Double.parseDouble(v[j]);
                }
                i+=1;
            }//while
        }catch(IOException ex){
            throw new RuntimeException(ex);
        }


        degK_k_g_eVecs = "degK_k_g_eVecs33";
        fileDegK_k_g_eVecs = new File(degK_k_g_eVecs);
        if(!fileDegK_k_g_eVecs.exists()){
            throw new RuntimeException("A file Does Not Exist!");
        }
        double[][] D;
        try{
            FileReader readerDegK_k_g_eVecs = new FileReader(degK_k_g_eVecs);
            BufferedReader buffDegK_k_g_eVecs = new BufferedReader(readerDegK_k_g_eVecs);
            String line;
            int nRows = 0;
            while((line = buffDegK_k_g_eVecs.readLine()) != null){
                nRows += 1;
            }
            int nCells = nC*nC*nC;
            int nColumns = 4+3*nBasis*(1+2*nBasis*3);
            D = new double[nRows][nColumns];
            int i = 0;
            readerDegK_k_g_eVecs = new FileReader(degK_k_g_eVecs);
            buffDegK_k_g_eVecs = new BufferedReader(readerDegK_k_g_eVecs);
            while((line = buffDegK_k_g_eVecs.readLine()) != null){
                line = line.trim();
                String[] v = line.split(" +");
                for (int j=0;j<nColumns;j++){
                    D[i][j] = Double.parseDouble(v[j]);
                }
                i+=1;
            }//while
        }catch(IOException ex){
            throw new RuntimeException(ex);
        }

        String c_modal_filename = "c_modal";
        File fileCModal = new File(c_modal_filename);

        if(!fileCModal.exists()){
            throw new RuntimeException("A file Does Not Exist!");
        }
        double[][] C_modal;
        try{
            FileReader readerCModal = new FileReader(c_modal_filename);
            BufferedReader buffCModal = new BufferedReader(readerCModal);
            String line;
            int nRows = 0;
            while((line = buffCModal.readLine()) != null){
                nRows += 1;
            }
            int nCells = nC*nC*nC;
            int nColumns = 2*3*nBasis*3*nBasis;
            C_modal = new double[nRows][nColumns];
            int i = 0;
            readerCModal = new FileReader(c_modal_filename);
            buffCModal = new BufferedReader(readerCModal);
            while((line = buffCModal.readLine()) != null){
                line = line.trim();
                String[] v = line.split(" +");
                for (int j=0;j<nColumns;j++){
                    C_modal[i][j] = Double.parseDouble(v[j]);
                }
                i+=1;
            }//while
        }catch(IOException ex){
            throw new RuntimeException(ex);
        }


        String c_modal11_filename = "c_modal11";
        File fileCModal11 = new File(c_modal11_filename);

        if(!fileCModal11.exists()){
            throw new RuntimeException("A file Does Not Exist!");
        }
        double[][] C_modal11;
        try{
            FileReader readerCModal = new FileReader(c_modal11_filename);
            BufferedReader buffCModal = new BufferedReader(readerCModal);
            String line;
            int nRows = 0;
            while((line = buffCModal.readLine()) != null){
                nRows += 1;
            }
            int nCells = nC*nC*nC;
            int nColumns = 2*3*nBasis*3*nBasis;
            C_modal11 = new double[nRows][nColumns];
            int i = 0;
            readerCModal = new FileReader(c_modal11_filename);
            buffCModal = new BufferedReader(readerCModal);
            while((line = buffCModal.readLine()) != null){
                line = line.trim();
                String[] v = line.split(" +");
                for (int j=0;j<nColumns;j++){
                    C_modal11[i][j] = Double.parseDouble(v[j]);
                }
                i+=1;
            }//while
        }catch(IOException ex){
            throw new RuntimeException(ex);
        }


        String c_modal22_filename = "c_modal22";
        File fileCModal22 = new File(c_modal22_filename);

        if(!fileCModal22.exists()){
            throw new RuntimeException("A file Does Not Exist!");
        }
        double[][] C_modal22;
        try{
            FileReader readerCModal = new FileReader(c_modal22_filename);
            BufferedReader buffCModal = new BufferedReader(readerCModal);
            String line;
            int nRows = 0;
            while((line = buffCModal.readLine()) != null){
                nRows += 1;
            }
            int nCells = nC*nC*nC;
            int nColumns = 2*3*nBasis*3*nBasis;
            C_modal22 = new double[nRows][nColumns];
            int i = 0;
            readerCModal = new FileReader(c_modal22_filename);
            buffCModal = new BufferedReader(readerCModal);
            while((line = buffCModal.readLine()) != null){
                line = line.trim();
                String[] v = line.split(" +");
                for (int j=0;j<nColumns;j++){
                    C_modal22[i][j] = Double.parseDouble(v[j]);
                }
                i+=1;
            }//while
        }catch(IOException ex){
            throw new RuntimeException(ex);
        }


        String c_modal33_filename = "c_modal33";
        File fileCModal33 = new File(c_modal33_filename);

        if(!fileCModal33.exists()){
            throw new RuntimeException("A file Does Not Exist!");
        }
        double[][] C_modal33;
        try{
            FileReader readerCModal = new FileReader(c_modal33_filename);
            BufferedReader buffCModal = new BufferedReader(readerCModal);
            String line;
            int nRows = 0;
            while((line = buffCModal.readLine()) != null){
                nRows += 1;
            }
            int nCells = nC*nC*nC;
            int nColumns = 2*3*nBasis*3*nBasis;
            C_modal33 = new double[nRows][nColumns];
            int i = 0;
            readerCModal = new FileReader(c_modal33_filename);
            buffCModal = new BufferedReader(readerCModal);
            while((line = buffCModal.readLine()) != null){
                line = line.trim();
                String[] v = line.split(" +");
                for (int j=0;j<nColumns;j++){
                    C_modal33[i][j] = Double.parseDouble(v[j]);
                }
                i+=1;
            }//while
        }catch(IOException ex){
            throw new RuntimeException(ex);
        }

        System.out.println("Running Lennard-Jones simulation");
        System.out.println(numAtoms+" atoms at density "+density+" and temperature "+temperature);
        System.out.println(" strain_x: " + strain_x);
        System.out.println(numSteps+" MC steps");
        final SimLJPropsNormalModes sim = new SimLJPropsNormalModes(Space.getInstance(3), numAtoms, density, temperature, rC, isLRC, strain_x);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE.setBox(sim.box);
        double ULat = meterPE.getDataAsScalar();
        MeterPressure  meterP = new MeterPressure(sim.space);
        meterP.setBox(sim.box);
        meterP.setTemperature(0);//ZERO means NO ideal gas component!
        meterP.setPotentialMaster(sim.potentialMaster);
        double PLat = meterP.getDataAsScalar();
        meterP.setTemperature(temperature);
        System.out.println(" Lattice Energy/N = " + (ULat/numAtoms) + " Lattice Pressure = " + PLat);


        double volume = sim.box.getBoundary().volume();

        MeterSolidPropsLJNormalModes meterUP = new MeterSolidPropsLJNormalModes(sim.space, new MeterPotentialEnergyFromIntegrator(sim.integrator), sim.potentialMaster, sim.coordinateDefinition, sim.cells0, temperature, dP, ULat, PLat,
                A, B, C, D, C_modal, C_modal11, C_modal22, C_modal33);
        A = null;
        B = null;
        C = null;
        D = null;
        C_modal = null;
        C_modal11 = null;
        C_modal22 = null;
        C_modal33 = null;

        //Initialization
        System.out.flush();
        long Ninit = numSteps/5;
        sim.initialize(Ninit);

        int numBlocks = 100;
        int interval = numAtoms;
        long blockSize = numSteps/(numBlocks*interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size "+blockSize+" interval "+interval);

        //U
        AccumulatorAverageFixed accumulatorPE = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPEPump = new DataPumpListener(meterPE, accumulatorPE, interval);
        sim.integrator.getEventManager().addListener(accumulatorPEPump);
        //P
        AccumulatorAverageFixed accumulatorP = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPPump = new DataPumpListener(meterP, accumulatorP, interval);
        sim.integrator.getEventManager().addListener(accumulatorPPump);

        //Mapped
//        AccumulatorAverageCovariance accumulatorUP = new AccumulatorAverageCovariance(blockSize);
        AccumulatorAverageFixed accumulatorUP = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorUPPump = new DataPumpListener(meterUP, accumulatorUP, interval);
        sim.integrator.getEventManager().addListener(accumulatorUPPump);

        final long startTime = System.currentTimeMillis();
/** RUN ... */
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        //U

        DataGroup dataPE = (DataGroup)accumulatorPE.getData();
        IData dataPEAvg = dataPE.getData(accumulatorPE.AVERAGE.index);
        IData dataPEErr = dataPE.getData(accumulatorPE.ERROR.index);
        IData dataPECorrelation = dataPE.getData(accumulatorPE.BLOCK_CORRELATION.index);
        double peAvg = dataPEAvg.getValue(0);
        double peErr = dataPEErr.getValue(0);
        double peCor = dataPECorrelation.getValue(0);
//        System.out.println(" Udirect/N = "+peAvg/numAtoms+"  Err = "+peErr/numAtoms+"  cor: "+peCor);
        //P
        DataGroup dataP = (DataGroup)accumulatorP.getData();
        IData dataPAvg = dataP.getData(accumulatorP.AVERAGE.index);
        IData dataPErr = dataP.getData(accumulatorP.ERROR.index);
        IData dataPCorrelation = dataP.getData(accumulatorP.BLOCK_CORRELATION.index);
        double pAvg = dataPAvg.getValue(0);
        double pErr = dataPErr.getValue(0);
        double pCor = dataPCorrelation.getValue(0);
        System.out.println(" Pdirect = "+pAvg+"  Err = "+pErr+"  cor: "+pCor);


        DataGroup dataUP = (DataGroup)accumulatorUP.getData();
        IData dataUPAvg = dataUP.getData(accumulatorUP.AVERAGE.index);
        IData dataUPErr = dataUP.getData(accumulatorUP.ERROR.index);
        IData dataUPCorr = dataUP.getData(accumulatorUP.BLOCK_CORRELATION.index);
//        IData dataUPCov = dataUP.getData(accumulatorUP.BLOCK_COVARIANCE.index);

        //P
        System.out.println(" Pconv " + dataUPAvg.getValue(0) +  "  "  + dataUPErr.getValue(0) +  "  " + dataUPCorr.getValue(0));
        System.out.println(" Phma  " + dataUPAvg.getValue(1) +  "  "  + dataUPErr.getValue(1) +  "  " + dataUPCorr.getValue(1));
        System.out.println(" Pnm   " + dataUPAvg.getValue(2) +  "  "  + dataUPErr.getValue(2) +  "  " + dataUPCorr.getValue(2));
        //P11
        System.out.println();
        System.out.println(" Pconv11 " + dataUPAvg.getValue(3) + "  " + dataUPErr.getValue(3) +  "  " + dataUPCorr.getValue(3));
        System.out.println(" Pconv22 " + dataUPAvg.getValue(4) + "  " + dataUPErr.getValue(4) +  "  " + dataUPCorr.getValue(4));
        System.out.println(" Pconv33 " + dataUPAvg.getValue(5) + "  " + dataUPErr.getValue(5) +  "  " + dataUPCorr.getValue(5));
        System.out.println();

        System.out.println(" Pnm11   " + dataUPAvg.getValue(6) + "  " + dataUPErr.getValue(6) + "   "+ dataUPCorr.getValue(6));
        System.out.println(" Pnm22   " + dataUPAvg.getValue(7) + "  " + dataUPErr.getValue(7) + "   "+ dataUPCorr.getValue(7));
        System.out.println(" Pnm33   " + dataUPAvg.getValue(8) + "  " + dataUPErr.getValue(8) + "   "+ dataUPCorr.getValue(8));
        System.out.println();
        System.out.println(" trace: 1/3(P11+P22+P33)");
        System.out.println(" Pconv_sum " + dataUPAvg.getValue(9) + "  " + dataUPErr.getValue(9) + "   "+ dataUPCorr.getValue(9));
        System.out.println(" Pnm_sum   " + dataUPAvg.getValue(10) + "  " + dataUPErr.getValue(10) + "   "+ dataUPCorr.getValue(10));

// Covariance
//        System.out.println("************************************************************************");
//        System.out.println(" Pvir11  " + dataUPAvg.getValue(0) + "  +/- " + dataUPErr.getValue(0) + "  cor: "+ dataUPCorr.getValue(0));
//        System.out.println(" Pvir22  " + dataUPAvg.getValue(1) + "  +/- " + dataUPErr.getValue(1) + "  cor: "+ dataUPCorr.getValue(1));
//        System.out.println(" Pnm11  " + dataUPAvg.getValue(2) + "  +/- " + dataUPErr.getValue(2) + "  cor: "+ dataUPCorr.getValue(2));
//        System.out.println(" Pnm22  " + dataUPAvg.getValue(3) + "  +/- " + dataUPErr.getValue(3) + "  cor: "+ dataUPCorr.getValue(3));
//        System.out.println("************************************************************************");
//
//        System.out.println(" conv1-conv2 " + dataUPCov.getValue(1)/Math.sqrt(dataUPCov.getValue(0)*dataUPCov.getValue(5)));
//        System.out.println(" nm1-nm2     " + dataUPCov.getValue(11)/Math.sqrt(dataUPCov.getValue(10)*dataUPCov.getValue(15)));

        long endTime = System.currentTimeMillis();
        System.out.println("time: " + (endTime - startTime)/1000.0);
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int nC = 5;
        public int nBasis = 4;
        public int numAtoms = nBasis*nC*nC*nC;
        public double density =  1.0;
        public double temperature  = 0.1;
        public boolean isLRC = false;
        public double dbP  = 1;
        public double rC = 3.0;
        public long numSteps = 1000000;
        public double strain_x = 0.0;
    }
}