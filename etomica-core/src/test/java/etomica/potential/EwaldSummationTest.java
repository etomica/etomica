package etomica.potential;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationResourceFile;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.models.water.P2WaterSPCE;
import etomica.models.water.SpeciesWater3P;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.units.Kelvin;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by Navneeth on 7/12/2017.
 */
public class EwaldSummationTest {

    EwaldSummation es;
    Box box;
    Simulation sim;

    public class ChargeAgentSourceSPCE implements AtomLeafAgentManager.AgentSource<EwaldSummation.MyCharge>{
        protected final EwaldSummation.MyCharge[] myCharge;
        public ChargeAgentSourceSPCE(SpeciesWater3P species){
            myCharge = new EwaldSummation.MyCharge[2];
            double chargeH = P2WaterSPCE.QH;
            double chargeO = P2WaterSPCE.QO;
            myCharge[species.getHydrogenType().getChildIndex()] = new EwaldSummation.MyCharge(chargeH);
            myCharge[species.getOxygenType().getChildIndex()] = new EwaldSummation.MyCharge(chargeO);
        }

        // *********************** set half(even # of particles ) as +ion, the other half -ion ***********************
        public EwaldSummation.MyCharge makeAgent(IAtom a, Box agentBox) {
            int index = a.getType().getChildIndex();
            return myCharge[index];
        }
        public void releaseAgent(EwaldSummation.MyCharge agent, IAtom atom, Box agentBox) {
            // Do nothing
        }

    }
    @Before
    public void setup(){
        Space space = Space.getInstance(3);
        box = new Box(space);
        SpeciesWater3P species = new SpeciesWater3P(space,false);
        ChargeAgentSourceSPCE agentSource = new ChargeAgentSourceSPCE(species);
        AtomLeafAgentManager<EwaldSummation.MyCharge> atomAgentManager = new AtomLeafAgentManager<EwaldSummation.MyCharge>(agentSource, box,EwaldSummation.MyCharge.class);
        double kcut = 5;
        double rCutRealES = 10;
        es = new EwaldSummation(box,atomAgentManager,space,kcut,rCutRealES);
        es.setAlpha(5.6/30);
        sim = new Simulation(space);

        sim.addSpecies(species);
        sim.addBox(box);
        box.setNMolecules(species,750);
        box.getBoundary().setBoxSize(new Vector3D(30,30,30));

        Configuration config = new ConfigurationResourceFile(
                String.format("etomica/potential/spce4.pos"),
                EwaldSummationTest.class
        );

        config.initializeCoordinates(box);



    }
    public static void main(String[] str) {
        EwaldSummationTest foo = new EwaldSummationTest();
        foo.setup();
        SimulationGraphic simg = new SimulationGraphic(foo.sim, foo.sim.getSpace(),null);
        simg.add(new DisplayBox(foo.sim, foo.box, foo.sim.getSpace(),null));
        simg.makeAndDisplayFrame();

    }

    @Test
    public void testUReal() throws Exception {
        double testval = es.uReal();
        double shouldbe = Kelvin.UNIT.toSim(-3.57226E+06);

        assertEquals(shouldbe,testval,1);

    }
    @Test
    public void testUFourier() throws Exception {
        double testval = es.uFourier();
        double shouldbe = Kelvin.UNIT.toSim(7.58785E+03);

        assertEquals(shouldbe,testval,0.01);

    }
    @Test
    public void testUSelf() throws Exception {
        double testval = es.uSelf();
        double shouldbe = Kelvin.UNIT.toSim(-1.42235E+07);

        assertEquals(shouldbe,testval,100);

    }

}
