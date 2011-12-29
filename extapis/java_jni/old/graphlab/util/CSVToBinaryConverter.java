package graphlab.util;
 
import java.io.*;


/**
 * Loads a comma-delimited file and makes a graph. Used by PageRank application.
 * @author akyrola
 *         Date: Oct 18, 2009
 */
public class CSVToBinaryConverter {

    public static void convertStochasticCSVToBinary(File f, File outfile) throws IOException {

        BufferedReader rd = new BufferedReader(new FileReader(f), 128000);
        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile), 128000));

        String ln;
        ln = rd.readLine();
        String[] tok = ln.split(",");

        assert(tok[0].equals("params"));
        int graphN = Integer.parseInt(tok[1]);

        dos.writeInt(graphN);



        /* Edges */
        while((ln = rd.readLine()) != null) {
            tok = ln.split(",");
            int i = Integer.parseInt(tok[0]);
            int j = Integer.parseInt(tok[1]);
            double weight = Double.parseDouble(tok[2]);
            dos.writeInt(i-1);
            dos.writeInt(j-1);
            dos.writeDouble(weight);
        }

        dos.writeInt(-1);
        dos.close();
    }

    public static void main(String[] args) throws IOException {
        File infile = new File(args[0]);
        File outfile = new File(args[1]);
        convertStochasticCSVToBinary(infile, outfile);
    }

}
