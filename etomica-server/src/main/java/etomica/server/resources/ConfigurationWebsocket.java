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
import java.util.concurrent.ScheduledFuture;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import static etomica.server.EtomicaServer.objectWriter;

@ServerEndpoint(
        value="/simulations/{id}/configuration",
        encoders = {ConfigurationWebsocket.ConfigurationUpdateEncoder.class}
)
@Metered
@Timed
public class ConfigurationWebsocket {
    private final SimulationStore simStore;
    private final ObjectMapper mapper;
    private final ScheduledThreadPoolExecutor executor;

    private final Logger log = LoggerFactory.getLogger(ConfigurationWebsocket.class);

    @Inject
    public ConfigurationWebsocket(SimulationStore store, ObjectMapper mapper, ScheduledThreadPoolExecutor executor) {
        this.simStore = store;
        this.mapper = mapper;
        this.executor = executor;
    }

    @OnOpen
    public void onOpen(final Session session, @PathParam("id") String id) {
        session.setMaxIdleTimeout(0);

        SimulationModel model = simStore.get(UUID.fromString(id));
        Simulation sim = model.getSimulation();
        SimulationWrapper wrapper = (SimulationWrapper) model.getWrapper(sim);

        Runnable sendConfigurationUpdate = () -> {
            if(sim.getController().isPaused() || !sim.getController().isActive()) {
                return;
            }


            sim.getController2().submitActionInterrupt(() -> {
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

        };

        ScheduledFuture<?> task = executor.scheduleWithFixedDelay(sendConfigurationUpdate, 0, 33, TimeUnit.MILLISECONDS);
        session.getUserProperties().put("task", task);
    }

    @OnClose
    public void onClose(Session session) {
        log.warn("Closing websocket");
        ((ScheduledFuture<?>) session.getUserProperties().get("task")).cancel(false);
    }

    @OnError
    public void onError(Session session, Throwable reason) {
        log.warn("Error in websocket", reason);
        ((ScheduledFuture<?>) session.getUserProperties().get("task")).cancel(false);
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
