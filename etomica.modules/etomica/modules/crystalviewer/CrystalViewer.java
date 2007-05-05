package etomica.modules.crystalviewer;

import java.awt.BorderLayout;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.JMenuItem;

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
import etomica.space.IVector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;


public class CrystalViewer {
    
    protected final Simulation sim;
    protected final JPanel panel;
    
    public CrystalViewer() {
        sim = new Simulation(Space3D.getInstance());
        sim.getDefaults().makeLJDefaults();
        center = sim.getSpace().makeVector();
        phase  = new Phase(sim);
        sim.addPhase(phase);

        species = new SpeciesSpheresMono(sim);
        sim.getSpeciesManager().addSpecies(species);
        
        panel = new JPanel(new BorderLayout());
        
        BasisMonatomic basisMonatomic = new BasisMonatomic(sim.getSpace());
        
        BravaisLattice[] lattices = new BravaisLattice[] {
                new LatticeCubicSimple(),
                new BravaisLatticeCrystal(new PrimitiveTetragonal(sim.getSpace()), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveHexagonal(sim.getSpace()), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveOrthorhombic(sim.getSpace()), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveMonoclinic(sim.getSpace()), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveTriclinic(sim.getSpace()), basisMonatomic),
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

        JPanel mainPanel = new JPanel();
        mainPanel.add(controlPanel);
        mainPanel.add(displayPhase.graphic());

        addMenu();
        panel.add(mainPanel);
    }

    private void addMenu() {
    	JMenuBar mBar = new JMenuBar();
    	JMenu fileMenu = new JMenu("File");
    	JMenuItem exitBtn = new JMenuItem("Exit");
    	exitBtn.addActionListener(new java.awt.event.ActionListener() {
    		public void actionPerformed(java.awt.event.ActionEvent ev) {
    			System.exit(0);
    		}
    	});
    	fileMenu.add(exitBtn);
    	JMenu helpMenu = new JMenu("Help");
    	JMenuItem aboutBtn = new JMenuItem("About Crystal Viewer");
    	aboutBtn.setEnabled(false);
    	helpMenu.add(aboutBtn);

    	mBar.add(fileMenu);
    	mBar.add(helpMenu);

    	panel.add(mBar, BorderLayout.NORTH);


    }

    public void update(BravaisLattice currentLattice) {
        latticePlane.setPrimitive(currentLattice.getPrimitive());
        displayPhase.setLabel(currentLattice.toString());
        clipPlaneEditor.update();
        displayPhase.repaint();
    }    
    
    public static class Applet extends javax.swing.JApplet {

	    public void init() {
            CrystalViewer viewer = new CrystalViewer();
		    getContentPane().add(viewer.panel);
	    }
    }
    
 
    public static void main(String[] args) {
        
        CrystalViewer viewer = new CrystalViewer();
        SimulationGraphic.makeAndDisplayFrame(viewer.panel, "Crystal Viewer");
    }
    
    protected SpeciesSpheresMono species;
    protected Phase phase;
    protected IVector center;
    protected LatticePlane latticePlane;
    protected ClipPlaneEditor clipPlaneEditor;
    protected DisplayPhase displayPhase;
    protected LatticeEditor latticeEditor;
        
}