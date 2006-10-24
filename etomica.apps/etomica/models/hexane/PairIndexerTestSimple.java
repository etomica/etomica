package etomica.models.hexane;

import java.util.ArrayList;
import java.util.HashMap;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomTypeSphere;
import etomica.atom.iterator.ApiInnerFixed;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.config.ConfigurationLattice;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorPhase;
import etomica.lattice.BravaisLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

public class PairIndexerTestSimple extends Simulation {

    public PairIndexerTestSimple(Space space, int latticeConfig){
        super(space);
      
        int numMolecules = -1;
        ConfigurationLattice config = null;
        defaults.makeLJDefaults();
        if (latticeConfig == SIMPLE_CUBIC) {
            numMolecules = 216;
            defaults.boxSize = 6;
            prim = new PrimitiveCubic(space);
            config = new ConfigurationLattice(new BravaisLattice(prim));
        }
        else if (latticeConfig == FCC) {
            numMolecules = 256;
            defaults.boxSize = 8;
            LatticeCubicFcc latticeFcc = new LatticeCubicFcc();
            config = new ConfigurationLattice(latticeFcc);
            prim = latticeFcc.getPrimitiveFcc();
        }
        defaults.atomSize = 0.5;
        
        config.setRememberingIndices(true);
        Species species = new SpeciesSpheresMono(this);
        phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numMolecules);
        BoundaryRectangularPeriodic bdry = new BoundaryRectangularPeriodic(this);
        
        phase.setBoundary(bdry);
        config.initializeCoordinates(phase);
        if (latticeConfig != SIMPLE_CUBIC) {
            ConfigurationLattice.MyLattice myLattice = (ConfigurationLattice.MyLattice)config.getLatticeMemento();
            Vector scaling = myLattice.latticeScaling;
            // assume isotropic scaling.  If it's not we're screwed anyway
            // because the FCC primitive has to be cubic
            prim.scaleSize(scaling.x(0));
        }
    }
    
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        
        PairIndexerTestSimple pit = new PairIndexerTestSimple(Space3D.getInstance(), FCC);
        
//      nan it will need to be changed*/
        PairIndexerMolecule pi = new PairIndexerMolecule(pit.phase, pit.prim);
      
        OutputFile printer = new OutputFile("Simple.txt");      
      
        AtomIterator inner = new AtomIteratorLeafAtoms(pit.phase);
        AtomIterator outer = new AtomIteratorLeafAtoms(pit.phase);
        ApiInnerFixed api = new ApiInnerFixed(outer, inner);
        api.reset();
        Vector work = pit.space.makeVector();

//      printer.println("Atom0:0 Atom0:1 Atom0:2 Atom0:# Atom1:0 Atom1:1 Atom1:2 Atom1:# Bin:");
      
        HashMap binHash = new HashMap();
        int maxBin = -1;
      
        while(api.hasNext()){
            AtomPair ap = api.nextPair();
            Atom atom0 = ap.atom0;
            Atom atom1 = ap.atom1;
            int bin = pi.getBin(atom0, atom1);

//        printer.println(atom0+" "+atom1+" "+dr+" "+bin);
//        printer.println(pi.getIndex(atom0)[0] +" "+ pi.getIndex(atom0)[1] + " " +
//                pi.getIndex(atom0)[2] + " " + atom0.node.getIndex() + " " +
//                pi.getIndex(atom1)[0] + " " + pi.getIndex(atom1)[1] + " " +
//                pi.getIndex(atom1)[2] + " " + atom1.node.getIndex() + " " +
//                pi.getBin(ap) +" " + atom0.getGlobalIndex());

          // determine the pair's dr
            Vector dr = pit.space.makeVector();
            dr.Ev1Mv2(((AtomLeaf)atom0).coord.position(), ((AtomLeaf)atom1).coord.position());
            if (bin < 0) {
                throw new RuntimeException("bin was negative, dr="+dr);
            }
            if (bin > maxBin) {
                maxBin = bin;
            }
            pit.phase.getBoundary().nearestImage(dr);
            if (bin != pi.getBin(atom1, atom0)) {
                pi.getBin(atom1, atom0);
                pi.getBin(atom0, atom1);
                System.out.println("asymmetry "+atom0+" "+atom1);
            }
          
            Integer iBin = new Integer(bin);
            // check if we've seen this bin before
            ArrayList drList = (ArrayList)binHash.get(iBin);
            if (drList == null) {
                // we haven't.  Add a new list for the bin containing our dr
                drList = new ArrayList();
                drList.add(dr);
                binHash.put(iBin, drList);
            }
            else {
                // we've seen it before.  Check to see if our dr is alread in 
                // the list
                boolean found = false;
                for (int i=0; i<drList.size(); i++) {
                    work.Ev1Mv2(dr, (Vector)drList.get(i));
                    if (work.squared() < 0.001) {
                        found = true;
                    }
                }
                if (!found) {
                    // it wasn't there, so add it
                    drList.add(dr);
                }
            }
        }
      
        // output a list of bins and all degenerate but apparently different vectors
        // for each bin.
        for (int i=0; i<maxBin; i++) {
            ArrayList drList = (ArrayList)binHash.get(new Integer(i));
            printer.println(i);
            if (drList != null) {
                for (int j=0; j<drList.size(); j++) {
                    printer.println(drList.get(j).toString());
                }
            }
        }
 
        printer.close();
 
        /*
        // uncomment to show the phase
        IntegratorPhase dummy = new IntegratorMC(pit.potentialMaster,1.0);
        dummy.setPhase(pit.phase);
        pit.getController().addAction(new ActivityIntegrate(dummy, false, false));
        ((AtomTypeSphere)pit.phase.firstAtom().type).setDiameter(1.0);
        SimulationGraphic simGraphic = new SimulationGraphic(pit);
        simGraphic.makeAndDisplayFrame();
        */
        
    }

    public Phase phase;
    public Primitive prim;
    public static final int SIMPLE_CUBIC = 0;
    public static final int FCC = 1;
}
