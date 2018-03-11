package etomica.parser.schroeder;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonSetter;
import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

@JsonIgnoreProperties(ignoreUnknown = true)
public class MDConfig {

    public int N;
    public int bondCount;
    public double size;
    public double gravity;
    public double dt;

    /**
     * format:
     * [ x, y, vx, vy, x2, y2, vx2, vy2, ... ]
     */
    @JsonSetter
    public void setData(double[] flatData) {
        data = new double[flatData.length / 4][];
        for (int i = 0; i < flatData.length; i+=4) {
            data[i / 4] = new double[]{
                    flatData[i],
                    flatData[i + 1],
                    flatData[i + 2],
                    flatData[i + 3]
            };
        }
    }

    public double[][] data;

    public int[] fixedList = new int[0];

    /**
     * every 2 indices is a bonded pair, [1a, 1b, 2a, 2b, ...]
     */
    @JsonSetter
    public void setBondList(int[] flatBondList) {
        if (flatBondList == null) {
            bondList = new int[0][];
        } else {
            bondList = new int[flatBondList.length / 2][];
            for (int i = 0; i < flatBondList.length; i+=2) {
                bondList[i / 2] = new int[]{
                        flatBondList[i],
                        flatBondList[i + 1]
                };
            }
        }
    }

    public int[][] bondList;

    public static void main(String[] args) throws IOException {
        ObjectMapper om = new ObjectMapper();
        om.configure(JsonParser.Feature.ALLOW_UNQUOTED_FIELD_NAMES, true);
        MDConfig config = om.readValue(new File("md.json"), MDConfig.class);
    }

}
