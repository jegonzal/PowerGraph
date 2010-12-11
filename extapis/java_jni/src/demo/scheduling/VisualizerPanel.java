package demo.scheduling;

import graphlab.Graph;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;

import demo.scheduling.ColorVertex;

/**
 * Shows a heatmap of vertices. Used in the demos.
 * @author akyrola
 *         Date: Jan 31, 2010
 */
public class VisualizerPanel extends JPanel {

    private int width;
    private int height;

    private BufferedImage img;
    private float recentTimeMax = 1.0f;
    private String mode = "Waiting...";

    public VisualizerPanel(int imgwidth, int imgheight, int width, int height) {
        this.width = width;
        this.height = height;
        img = new BufferedImage(imgwidth, imgheight, BufferedImage.TYPE_INT_BGR);
        img.createGraphics().setColor(Color.RED);
    }



    public Dimension getPreferredSize() {
        return new Dimension(width, height);
    }

    public float getRecentTimeMax() {
        return recentTimeMax;
    }

    public void setRecentTimeMax(float recentTimeMax) {
        this.recentTimeMax = recentTimeMax;
    }

    public String getMode() {
        return mode;
    }

    public void setMode(String mode) {
        this.mode = mode;
    }

    public void updateData(Graph graph) {
        int[] scaledvals = new int[graph.getNumOfVertices()];
        for(int vid=0; vid<graph.getNumOfVertices(); vid++) {
            scaledvals[vid] = ((ColorVertex) graph.getVertex(vid)).getValue();
        }
        /* Update image  */
        int[] imgdata = ((DataBufferInt) img.getRaster().getDataBuffer()).getData();
        try {
            for(int i = 0; i<scaledvals.length; i++) {
                imgdata[i] = colors[(255-scaledvals[i])%256];
            }
        } catch (ArrayIndexOutOfBoundsException aio) {
            System.err.println("Array out of bounds:  mode=" + mode);
            throw aio;
        }
        this.repaint();
    }

    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        g.drawImage(img, 0, 0, this.getWidth(), this.getHeight()-50, null);
        //g.drawString("Mode: " + mode + "[" + info + "]", 10, this.getHeight()-30);
    }


    // Heatmap colors, color map from R-software: heat.colors(256)

    static int colors[] = { 0xFF0000FF, 0xFF0100FF, 0xFF0300FF, 0xFF0400FF, 0xFF0500FF, 0xFF0700FF,
            0xFF0800FF, 0xFF0900FF, 0xFF0B00FF, 0xFF0C00FF, 0xFF0D00FF, 0xFF0F00FF,
            0xFF1000FF, 0xFF1100FF, 0xFF1300FF, 0xFF1400FF, 0xFF1500FF, 0xFF1700FF,
            0xFF1800FF, 0xFF1900FF, 0xFF1B00FF, 0xFF1C00FF, 0xFF1D00FF, 0xFF1F00FF,
            0xFF2000FF, 0xFF2100FF, 0xFF2300FF, 0xFF2400FF, 0xFF2500FF, 0xFF2700FF,
            0xFF2800FF, 0xFF2900FF, 0xFF2B00FF, 0xFF2C00FF, 0xFF2D00FF, 0xFF2F00FF,
            0xFF3000FF, 0xFF3100FF, 0xFF3300FF, 0xFF3400FF, 0xFF3500FF, 0xFF3700FF,
            0xFF3800FF, 0xFF3900FF, 0xFF3B00FF, 0xFF3C00FF, 0xFF3D00FF, 0xFF3F00FF,
            0xFF4000FF, 0xFF4100FF, 0xFF4300FF, 0xFF4400FF, 0xFF4500FF, 0xFF4700FF,
            0xFF4800FF, 0xFF4900FF, 0xFF4B00FF, 0xFF4C00FF, 0xFF4D00FF, 0xFF4F00FF,
            0xFF5000FF, 0xFF5100FF, 0xFF5300FF, 0xFF5400FF, 0xFF5500FF, 0xFF5700FF,
            0xFF5800FF, 0xFF5900FF, 0xFF5B00FF, 0xFF5C00FF, 0xFF5D00FF, 0xFF5F00FF,
            0xFF6000FF, 0xFF6100FF, 0xFF6300FF, 0xFF6400FF, 0xFF6500FF, 0xFF6700FF,
            0xFF6800FF, 0xFF6900FF, 0xFF6B00FF, 0xFF6C00FF, 0xFF6D00FF, 0xFF6F00FF,
            0xFF7000FF, 0xFF7100FF, 0xFF7300FF, 0xFF7400FF, 0xFF7500FF, 0xFF7700FF,
            0xFF7800FF, 0xFF7900FF, 0xFF7B00FF, 0xFF7C00FF, 0xFF7D00FF, 0xFF7F00FF,
            0xFF8000FF, 0xFF8200FF, 0xFF8300FF, 0xFF8400FF, 0xFF8600FF, 0xFF8700FF,
            0xFF8800FF, 0xFF8A00FF, 0xFF8B00FF, 0xFF8C00FF, 0xFF8E00FF, 0xFF8F00FF,
            0xFF9000FF, 0xFF9200FF, 0xFF9300FF, 0xFF9400FF, 0xFF9600FF, 0xFF9700FF,
            0xFF9800FF, 0xFF9A00FF, 0xFF9B00FF, 0xFF9C00FF, 0xFF9E00FF, 0xFF9F00FF,
            0xFFA000FF, 0xFFA200FF, 0xFFA300FF, 0xFFA400FF, 0xFFA600FF, 0xFFA700FF,
            0xFFA800FF, 0xFFAA00FF, 0xFFAB00FF, 0xFFAC00FF, 0xFFAE00FF, 0xFFAF00FF,
            0xFFB000FF, 0xFFB200FF, 0xFFB300FF, 0xFFB400FF, 0xFFB600FF, 0xFFB700FF,
            0xFFB800FF, 0xFFBA00FF, 0xFFBB00FF, 0xFFBC00FF, 0xFFBE00FF, 0xFFBF00FF,
            0xFFC000FF, 0xFFC200FF, 0xFFC300FF, 0xFFC400FF, 0xFFC600FF, 0xFFC700FF,
            0xFFC800FF, 0xFFCA00FF, 0xFFCB00FF, 0xFFCC00FF, 0xFFCE00FF, 0xFFCF00FF,
            0xFFD000FF, 0xFFD200FF, 0xFFD300FF, 0xFFD400FF, 0xFFD600FF, 0xFFD700FF,
            0xFFD800FF, 0xFFDA00FF, 0xFFDB00FF, 0xFFDC00FF, 0xFFDE00FF, 0xFFDF00FF,
            0xFFE000FF, 0xFFE200FF, 0xFFE300FF, 0xFFE400FF, 0xFFE600FF, 0xFFE700FF,
            0xFFE800FF, 0xFFEA00FF, 0xFFEB00FF, 0xFFEC00FF, 0xFFEE00FF, 0xFFEF00FF,
            0xFFF000FF, 0xFFF200FF, 0xFFF300FF, 0xFFF400FF, 0xFFF600FF, 0xFFF700FF,
            0xFFF800FF, 0xFFFA00FF, 0xFFFB00FF, 0xFFFC00FF, 0xFFFE00FF, 0xFFFF00FF,
            0xFFFF02FF, 0xFFFF06FF, 0xFFFF0AFF, 0xFFFF0EFF, 0xFFFF12FF, 0xFFFF16FF,
            0xFFFF1AFF, 0xFFFF1EFF, 0xFFFF22FF, 0xFFFF26FF, 0xFFFF2AFF, 0xFFFF2EFF,
            0xFFFF32FF, 0xFFFF36FF, 0xFFFF3AFF, 0xFFFF3EFF, 0xFFFF42FF, 0xFFFF46FF,
            0xFFFF4AFF, 0xFFFF4EFF, 0xFFFF52FF, 0xFFFF56FF, 0xFFFF5AFF, 0xFFFF5EFF,
            0xFFFF62FF, 0xFFFF66FF, 0xFFFF6AFF, 0xFFFF6EFF, 0xFFFF72FF, 0xFFFF76FF,
            0xFFFF7AFF, 0xFFFF7EFF, 0xFFFF81FF, 0xFFFF85FF, 0xFFFF89FF, 0xFFFF8DFF,
            0xFFFF91FF, 0xFFFF95FF, 0xFFFF99FF, 0xFFFF9DFF, 0xFFFFA1FF, 0xFFFFA5FF,
            0xFFFFA9FF, 0xFFFFADFF, 0xFFFFB1FF, 0xFFFFB5FF, 0xFFFFB9FF, 0xFFFFBDFF,
            0xFFFFC1FF, 0xFFFFC5FF, 0xFFFFC9FF, 0xFFFFCDFF, 0xFFFFD1FF, 0xFFFFD5FF,
            0xFFFFD9FF, 0xFFFFDDFF, 0xFFFFE1FF, 0xFFFFE5FF, 0xFFFFE9FF, 0xFFFFEDFF,
            0xFFFFF1FF, 0xFFFFF5FF, 0xFFFFF9FF, 0xFFFFFDFF
    };

}
