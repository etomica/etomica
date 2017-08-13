package etomica.server.resources;

import com.codahale.metrics.annotation.Metered;
import com.codahale.metrics.annotation.Timed;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.meta.SimulationModel;
import etomica.simulation.Simulation;

import javax.websocket.*;
import javax.websocket.server.PathParam;
import javax.websocket.server.ServerEndpoint;
import java.util.Map;
import java.util.Timer;
import java.util.TimerTask;
import java.util.UUID;

@ServerEndpoint(
        value="/simulations/{id}/configuration",
        encoders = {ConfigurationStreamResource.ConfigurationEncoder.class}
)
@Metered
@Timed
public class ConfigurationStreamResource {


    @OnOpen
    public void onOpen(final Session session, @PathParam("id") String id) {
        Map<UUID, SimulationModel> simStore = (Map<UUID, SimulationModel>) session.getUserProperties().get("simStore");
        Timer timer = (Timer) session.getUserProperties().get("timer");
        SimulationModel model = simStore.get(UUID.fromString(id));

        timer.schedule(new ConfigurationTimerTask(session, model.getSimulation()), 0, 33);

    }

    private static class ConfigurationTimerTask extends TimerTask {
        private final Session session;
        private final Simulation simulation;

        private ConfigurationTimerTask(Session session, Simulation simulation) {
            this.session = session;
            this.simulation = simulation;
        }

        @Override
        public void run() {
            if(simulation.getController().isPaused() || !simulation.getController().isActive()) {
                return;
            }


            // an array of vector coordinates for each box
            double[][][] boxes = new double[simulation.getBoxCount()][][];
            simulation.getController().doActionNow(() -> {
                for (int i = 0; i < simulation.getBoxCount(); i++) {
                    Box box = simulation.getBox(i);
                    IAtomList atomList = box.getLeafList();
                    boxes[i] = new double[atomList.getAtomCount()][simulation.getSpace().D()];

                    for (int j = 0; j < atomList.getAtomCount(); j++) {
                        atomList.getAtom(j).getPosition().assignTo(boxes[i][j]);
                    }
                }
            });

            session.getAsyncRemote().sendObject(boxes);
        }
    }

    public static class ConfigurationEncoder implements Encoder.Text<double[][][]> {
        // TODO: inject this from application, should be ok for now
        private static final ObjectMapper mapper = new ObjectMapper();

        @Override
        public String encode(double[][][] object) throws EncodeException {
            try {
                return mapper.writeValueAsString(object);
            } catch (JsonProcessingException e) {
                throw new EncodeException(object, "Json error", e);
            }
        }

        @Override
        public void init(EndpointConfig config) {

        }

        @Override
        public void destroy() {

        }
    }
}
