package org.example;

import org.apache.commons.imaging.*;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

public class ImageUtils {
    public static double[][] loadImage(String path) throws Exception {
        BufferedImage image = Imaging.getBufferedImage(new File(path));
        // Convert to grayscale and resize if necessary
        return convertToGrayscale(image);
    }

    public static void writeImage(BufferedImage image, String filePath) throws IOException {
        ImageIO.write(image, "jpg", new File(filePath));
    }

    private static double[][] convertToGrayscale(BufferedImage image) {
        int width = image.getWidth();
        int height = image.getHeight();
        double[][] grayscale = new double[width][height];

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                Color color = new Color(image.getRGB(i, j));
                int red = color.getRed();
                int green = color.getGreen();
                int blue = color.getBlue();
                // Convert to grayscale value
                grayscale[i][j] = 0.299 * red + 0.587 * green + 0.114 * blue;
                // Normalize the pixel values
                grayscale[i][j] -= 128;
            }
        }
        return grayscale;
    }

    public static double calculateMSE(BufferedImage original, BufferedImage reconstructed) {
        double mse = 0;
        for (int y = 0; y < original.getHeight(); y++) {
            for (int x = 0; x < original.getWidth(); x++) {
                int originalPixel = new Color(original.getRGB(x, y)).getRed();
                int reconstructedPixel = new Color(reconstructed.getRGB(x, y)).getRed();
                mse += Math.pow(originalPixel - reconstructedPixel, 2);
            }
        }
        mse /= (original.getWidth() * original.getHeight());
        return mse;
    }

    public static double calculatePSNR(double mse, double maxPixelValue) {
        return 20 * Math.log10(maxPixelValue / Math.sqrt(mse));
    }

}
