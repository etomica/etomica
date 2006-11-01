package etomica.modules.crystalviewer;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.BravaisLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeCubicBcc;
import etomica.lattice.LatticeCubicDiamond;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.lattice.LatticeHcp;
import etomica.lattice.LatticePlane;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.lattice.crystal.PrimitiveMonoclinic;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.lattice.crystal.PrimitiveTriclinic;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;


public class CrystalViewer {
    
    protected final Simulation sim;
    protected final JPanel panel;
    
    public CrystalViewer() {
        sim = new Simulation(Space3D.getInstance());
        sim.getDefaults().makeLJDefaults();
        center = sim.space.makeVector();
        phase  = new Phase(sim);

        species = new SpeciesSpheresMono(sim);
        
        panel = new JPanel();
        
        double a = 1.0;//default lattice constant
        BasisMonatomic basisMonatomic = new BasisMonatomic(sim.space);
        
        BravaisLattice[] lattices = new BravaisLattice[] {
                new LatticeCubicSimple(),
                new BravaisLatticeCrystal(new PrimitiveTetragonal(sim.space), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveHexagonal(sim.space), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveOrthorhombic(sim.space), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveMonoclinic(sim.space), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveTriclinic(sim.space), basisMonatomic),
                new LatticeCubicFcc(),
                new LatticeCubicBcc(),
                new LatticeHcp(),
                new LatticeCubicDiamond()
            };
        String[] latticeNames = new String[]{
                "Simple Cubic", "Tetragonal", "Hexagonal", "Orthorhombic", "Monoclinic", "Triclinic", "FCC", "BCC", "HCP", "Diamond"};

        displayPhase = new DisplayPhase(phase);
        displayPhase.setPixelUnit(new Pixel(20));

        // we pass these to make LatticePlane happy.  they'll get whacked by update() later
        latticePlane = new LatticePlane(lattices[0].getPrimitive(), new int[] {1,0,0});
        
        clipPlaneEditor = new ClipPlaneEditor(latticePlane, displayPhase);
        
        latticeEditor = new LatticeEditor(this, lattices, latticeNames);
        
        JTabbedPane controlTabs = new javax.swing.JTabbedPane();
        controlTabs.add("Crystal", latticeEditor.getPanel());
        controlTabs.add("Plane", clipPlaneEditor.getPanel());
        JPanel controlPanel = new JPanel();
        controlPanel.add(controlTabs);
        panel.add(controlPanel);
        panel.add(displayPhase.graphic());
    }
    
    public void update(BravaisLattice currentLattice) {
        latticePlane.setPrimitive(currentLattice.getPrimitive());
        displayPhase.setLabel(currentLattice.toString());
        clipPlaneEditor.update();
    }    
    
    public static class Applet extends javax.swing.JApplet {

	    public void init() {
            CrystalViewer viewer = new CrystalViewer();
		    getContentPane().add(viewer.panel);
	    }
    }
    
 
    public static void main(String[] args) {
        
        CrystalViewer viewer = new CrystalViewer();
        SimulationGraphic.makeAndDisplayFrame(viewer.panel);
    }
    
    protected SpeciesSpheresMono species;
    protected Phase phase;
    protected Vector center;
    protected LatticePlane latticePlane;
    protected ClipPlaneEditor clipPlaneEditor;
    protected DisplayPhase displayPhase;
    protected LatticeEditor latticeEditor;
        
}