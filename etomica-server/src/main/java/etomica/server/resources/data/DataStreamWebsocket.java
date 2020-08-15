package etomica.server.resources.data;

import com.codahale.metrics.annotation.Metered;
import com.codahale.metrics.annotation.Timed;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.data.*;
import etomica.meta.SimulationModel;
import etomica.server.dao.DataStreamStore;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.DataAndInfo;
import etomica.simulation.Simulation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.inject.Inject;
import javax.websocket.*;
import javax.websocket.server.PathParam;
import javax.websocket.server.ServerEndpoint;
import java.util.Timer;
import java.util.TimerTask;
import java.util.UUID;
import java.util.concurrent.ScheduledFuture;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

@ServerEndpoint(
        value="/simulations/{simId}/data/{dataId}",
        encoders = {DataStreamWebsocket.DataAndInfoEncoder.class}
)
@Metered
@Timed
public class DataStreamWebsocket {
    private final SimulationStore simStore;
    private final DataStreamStore dataStore;
    private final ObjectMapper mapper;
    private final ScheduledThreadPoolExecutor executor;

    private final Logger log = LoggerFactory.getLogger(DataStreamWebsocket.class);

    @Inject
    public DataStreamWebsocket(SimulationStore simStore, DataStreamStore dataStore, ObjectMapper mapper, ScheduledThreadPoolExecutor executor) {
        this.simStore = simStore;
        this.dataStore = dataStore;
        this.mapper = mapper;
        this.executor = executor;
    }

    @OnOpen
    public void onOpen(final Session session, @PathParam("simId") String simId, @PathParam("dataId") String dataId) {
        session.setMaxIdleTimeout(0);

        session.getUserProperties().put("mapper", mapper);
        SimulationModel model = simStore.get(UUID.fromString(simId));
        Simulation sim = model.getSimulation();
        DataStreamStore.DataPlumbing dataPlumbing = dataStore.get(UUID.fromString(dataId));
        DataDump dump = dataPlumbing.getDump();
        final DataAndInfo dataAndInfo = new DataAndInfo();

        Runnable sendData = () -> {
            if(sim.getController().isPaused() || !sim.getController().isActive()) {
                return;
            }

            sim.getController2().submitActionInterrupt(() -> {
                IData data = dump.getData();
                dataAndInfo.setData(dump.getData());
                dataAndInfo.setDataInfo(dump.getDataInfo());
                if(data != null) {
                    session.getAsyncRemote().sendObject(dataAndInfo);
                }
            });
        };

        ScheduledFuture<?> task = executor.scheduleWithFixedDelay(sendData, 0, 333, TimeUnit.MILLISECONDS);
        session.getUserProperties().put("task", task);
        // add on construction
//        model.getSimulation().getIntegrator().getEventManager().addListener(dataPlumbing.getPump());
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

    public static class DataAndInfoEncoder implements Encoder.Text<DataAndInfo> {
        private EndpointConfig config;

        @Override
        public String encode(DataAndInfo object) throws EncodeException {
            ObjectMapper mapper = (ObjectMapper) config.getUserProperties().get("mapper");
            try {
                return mapper.writeValueAsString(object);
            } catch (JsonProcessingException e) {
                throw new EncodeException(object, "JSON error", e);
            }
        }

        @Override
        public void init(EndpointConfig config) {
            this.config = config;
        }

        @Override
        public void destroy() {

        }
    }

}
