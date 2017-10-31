package etomica.server.resources;

import com.codahale.metrics.annotation.Metered;
import com.codahale.metrics.annotation.Timed;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectWriter;
import etomica.meta.SimulationModel;
import etomica.meta.wrappers.SimulationWrapper;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.ConfigurationUpdate;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.inject.Inject;
import javax.websocket.*;
import javax.websocket.server.PathParam;
import javax.websocket.server.ServerEndpoint;
import java.io.IOException;
import java.io.Writer;
import java.util.Timer;
import java.util.TimerTask;
import java.util.UUID;

import static etomica.server.EtomicaServer.objectWriter;

@ServerEndpoint(
        value="/simulations/{id}/configuration",
        encoders = {ConfigurationWebsocket.ConfigurationUpdateEncoder.class}
)
@Metered
@Timed
public class ConfigurationWebsocket {
    private final SimulationStore simStore;
    private final Timer timer;
    private final ObjectMapper mapper;

    private final Logger log = LoggerFactory.getLogger(ConfigurationWebsocket.class);

    @Inject
    public ConfigurationWebsocket(SimulationStore store, Timer timer, ObjectMapper mapper) {
        this.simStore = store;
        this.timer = timer;
        this.mapper = mapper;
    }

    @OnOpen
    public void onOpen(final Session session, @PathParam("id") String id) {
        SimulationModel model = simStore.get(UUID.fromString(id));

        timer.schedule(new ConfigurationTimerTask(session, model), 0, 33);
//        for (int i = 0; i < model.getSimulation().getBoxCount(); i++) {
//            model.getSimulation().getBox(i).getBoundary().getEventManager().addListener(
//                    new WSBoundaryListener(session, model, mapper)
//            );
//        }
    }

    @OnClose
    public void onClose(Session session) {
        log.warn("Closing websocket");
    }

    @OnError
    public void onError(Session session, Throwable reason) {
        log.warn("Error in websocket", reason);
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
                Boundary[] boundaries = new Boundary[sim.getBoxCount()];
                for (int i = 0; i < sim.getBoxCount(); i++) {
                    boundaries[i] = sim.getBox(i).getBoundary();
                }

                ConfigurationUpdate update = new ConfigurationUpdate(
                        wrapper.getAllCoordinates(),
                        boundaries
                );
                session.getAsyncRemote().sendObject(update);
            });
        }
    }

    public static class ConfigurationUpdateEncoder implements Encoder.TextStream<ConfigurationUpdate> {
        private static ObjectWriter objectWriter = objectWriter();

        @Override
        public void encode(ConfigurationUpdate object, Writer writer) throws EncodeException, IOException {
            objectWriter.writeValue(writer, object);

        }

        @Override
        public void init(EndpointConfig config) {

        }

        @Override
        public void destroy() {

        }
    }
}
