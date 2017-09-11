package etomica.server.resources;

import com.codahale.metrics.annotation.Metered;
import com.codahale.metrics.annotation.Timed;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.meta.SimulationModel;
import etomica.meta.wrappers.SimulationWrapper;
import etomica.server.dao.SimulationStore;
import etomica.simulation.Simulation;

import javax.inject.Inject;
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
    private final SimulationStore simStore;
    private final Timer timer;

    @Inject
    public ConfigurationStreamResource(SimulationStore store, Timer timer) {
        this.simStore = store;
        this.timer = timer;
    }

    @OnOpen
    public void onOpen(final Session session, @PathParam("id") String id) {
        SimulationModel model = simStore.get(UUID.fromString(id));

        timer.schedule(new ConfigurationTimerTask(session, model), 0, 33);

    }

    private static class ConfigurationTimerTask extends TimerTask {
        private final Session session;

        private final SimulationWrapper wrapper;
        private final Simulation sim;

        private ConfigurationTimerTask(Session session, SimulationModel model) {
            this.session = session;
            this.wrapper = (SimulationWrapper) model.getWrapper(model.getSimulation());
            this.sim = model.getSimulation();

        }

        @Override
        public void run() {
            if(sim.getController().isPaused() || !sim.getController().isActive()) {
                return;
            }


            sim.getController().doActionNow(() -> {
                session.getAsyncRemote().sendObject(wrapper.getAllCoordinates());
            });
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
