package etomica.server.dao;

import etomica.data.DataDump;
import etomica.data.DataPumpListener;

import javax.inject.Singleton;
import java.util.UUID;
import java.util.concurrent.ConcurrentHashMap;

@Singleton
public class DataStreamStore extends ConcurrentHashMap<UUID, DataStreamStore.DataPlumbing> {
    @Override
    public DataPlumbing get(Object o) {
        if (o.equals(UUID.fromString("11111111-1111-1111-1111-111111111111"))) {
            return (DataPlumbing) super.values().toArray()[0];
        }
        return super.get(o);
    }

    public static class DataPlumbing {
        private final DataPumpListener pump;
        private final DataDump dump;
        private final UUID simId;

        public DataPlumbing(DataPumpListener pump, DataDump dump, UUID simId) {
            this.pump = pump;
            this.dump = dump;
            this.simId = simId;
        }

        public DataPumpListener getPump() {
            return pump;
        }

        public DataDump getDump() {
            return dump;
        }

        public UUID getSimId() {
            return simId;
        }
    }
}

