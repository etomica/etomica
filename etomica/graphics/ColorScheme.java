package etomica.graphics;
import java.awt.Color;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.atom.IAtom;
import etomica.integrator.IntegratorHard;
import etomica.modifier.ModifierGeneral;
import etomica.nbr.list.PotentialMasterList;
import etomica.util.Default;

/**
 * Class that defines the algorithm used to determine atoms colors when drawn to DisplayPhase.
 * The atomColor method is called just before the atom is drawn to set the graphics color.
 *
 * @author David Kofke
 */
 
public abstract class ColorScheme implements java.io.Serializable {

    protected Color defaultColor;
    
    public ColorScheme() {
        this(DEFAULT_ATOM_COLOR);
    }
    public ColorScheme(Color color) {
        defaultColor = color;
    }

    public abstract Color getAtomColor(IAtom a);
    
    public final void setDefaultColor(Color c) {defaultColor = c;}
    public final Color getDefaultColor() {return defaultColor;}

    public static Color DEFAULT_ATOM_COLOR = Color.red;
    
    /**
     * Colors all atoms with baseColor.
     */
    public static class Simple extends ColorScheme {
        private static final long serialVersionUID = 1L;
        public Simple() {super();}
        public Simple(java.awt.Color color) {super(color);}
        public Color getAtomColor(IAtom a) {return defaultColor;}
    }
    
    //for color scheme demo
    public static void main(String args[]) {
      Default defaults = new Default();
      defaults.doSleep = false;
      defaults.ignoreOverlap = true;
      final etomica.simulation.prototypes.HSMD3D sim = new etomica.simulation.prototypes.HSMD3D();
      final SimulationGraphic simGraphic = new SimulationGraphic(sim);
      DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
      nSelector.setResetAction(new SimulationRestart(sim));
      nSelector.setSpeciesAgent(sim.phase.getAgent(sim.species));
      simGraphic.add(nSelector);
      
      final ColorSchemeByType ct = new ColorSchemeByType();
      final ColorSchemeTemperature ctemp = new ColorSchemeTemperature(0,5);
      final ColorSchemeColliders ccld =
        new ColorSchemeColliders((IntegratorHard)sim.getIntegratorList().getFirst());
      final ColorSchemeNeighbor nghb = new ColorSchemeNeighbor(sim, (PotentialMasterList)sim.potentialMaster, sim.phase);
      nghb.setAtom(sim.phase.getSpeciesMaster().getLeafList().getAtom(0));
      final ColorSchemeRandom rand = new ColorSchemeRandom(sim.phase, sim.getRandom());
      final ColorSchemeCell cell = new ColorSchemeCell((PotentialMasterList)sim.potentialMaster,sim.getRandom(),sim.phase);
      cell.setLattice(((PotentialMasterList)sim.potentialMaster).getNbrCellManager(sim.phase).getLattice());
      
      Action act = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          ct.setColor(sim.species.getMoleculeType(), 
              new java.awt.Color(sim.getRandom().nextInt(256),
                  sim.getRandom().nextInt(256),
                  sim.getRandom().nextInt(256))
              );
          dp.setColorScheme(ct);
        }
      };
      Action act2 = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(ctemp);
        }
      };
      Action act3 = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(ccld);
        }
      };
      Action act4 = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(nghb);
        }
      };
      Action act5 = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(rand);
        }
      };
      Action act6 = new Action() {
        public void actionPerformed() {
          ct.setColor(sim.species.getMoleculeType(), java.awt.Color.red);
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(ct);
        }
      };
      Action act7 = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(cell);
        }
      };

      DeviceButton colorer = new DeviceButton(sim.getController(),act);
      colorer.setLabel("Global Random");
      DeviceButton tempcolorer = new DeviceButton(sim.getController(),act2);
      tempcolorer.setLabel("By Temperature (0-5)");
      DeviceButton colliders = new DeviceButton(sim.getController(),act3);
      colliders.setLabel("Colliders");
      DeviceButton neighbors = new DeviceButton(sim.getController(),act4);
      neighbors.setLabel("Neighbors");
      DeviceButton randomcol = new DeviceButton(sim.getController(),act5);
      randomcol.setLabel("Unique Random");
      DeviceButton def = new DeviceButton(sim.getController(),act6);
      def.setLabel("Default Red");
      DeviceButton cellbtn = new DeviceButton(sim.getController(),act7);
      cellbtn.setLabel("Cell");
      
      DeviceSlider slabslide = new DeviceSlider(sim.getController());
      slabslide.setMinimum(0);
      slabslide.setMaximum(100);
      final DisplayPhaseCanvasG3DSys dpg3d = (DisplayPhaseCanvasG3DSys)simGraphic.getDisplayPhase(sim.phase).canvas;
      ModifierGeneral mg = new ModifierGeneral(dpg3d,"slab");
      slabslide.setModifier(mg);
      slabslide.setName("Slab");
      
      DeviceSlider depthslide = new DeviceSlider(sim.getController());
      depthslide.setMinimum(0);
      depthslide.setMaximum(100);
      ModifierGeneral mg2 = new ModifierGeneral(dpg3d,"depth");
      depthslide.setModifier(mg2);
      depthslide.setName("Depth");
      
      simGraphic.add(def);
      simGraphic.add(colorer);
      simGraphic.add(tempcolorer);
      simGraphic.add(colliders);
      simGraphic.add(neighbors);
      simGraphic.add(randomcol);
      simGraphic.add(cellbtn);
      simGraphic.add(slabslide);
      simGraphic.add(depthslide);
      
      
//      d.get
//      Action setslab = new Action() {
//
//        public void actionPerformed() {
//          //d.setS
//        }
//      };
      
      simGraphic.makeAndDisplayFrame();
      
      //ColorSchemeTemperature colorScheme = new ColorSchemeTemperature(0.0f,5.0f);

      ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simGraphic.displayList().getFirst()).getColorScheme());
      colorScheme.setColor(sim.species.getMoleculeType(), java.awt.Color.red);

      DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
      dp.setColorScheme(colorScheme);
      simGraphic.panel().setBackground(java.awt.Color.yellow);        
    }
}
