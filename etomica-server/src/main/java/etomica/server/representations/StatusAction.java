package etomica.server.representations;

import com.fasterxml.jackson.annotation.JsonProperty;

public class StatusAction {
    private Status status;

    public Status getStatus() {
        return status;
    }

    public void setStatus(Status status) {
        this.status = status;
    }

    public enum Status {
        @JsonProperty("start")
        START,

        @JsonProperty("pause")
        PAUSE,

        @JsonProperty("reset")
        RESET

    }
}

