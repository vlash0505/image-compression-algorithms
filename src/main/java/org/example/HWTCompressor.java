package org.example;

import org.apache.commons.imaging.ImageReadException;
import org.apache.commons.imaging.Imaging;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/*
Implementation of the Discrete Wavelet Transform (DWT).

Step 1:
Apply Haar Transform for rows and columns.

Step 2:
Compress the image data.

Step 3:
Apply Inverse Haar Transform for rows and columns.

Step 4:
Reconstruct the image.
*/

public class HWTCompressor {

    public static void main(String[] args) {
        try {
            String inputImagePath = "rhino.jpg";
            String outputImagePath = "result_hwt.jpg";
            double threshold = 2; // Set the threshold for compression

            BufferedImage image = readImage(inputImagePath);
            double[][] imageData = convertTo2DUsingGetRGB(image);

            // Apply Haar Transform
            applyHaarTransform(imageData, true);
            applyHaarTransform(imageData, false);

            // Compress the image data
            compressWaveletCoefficients(imageData, threshold);

            // Apply Inverse Haar Transform
            applyInverseHaarTransform(imageData, false);
            applyInverseHaarTransform(imageData, true);

            // Reconstruct and write the image
            BufferedImage outputImage = reconstructImage(imageData);
            writeImage(outputImage, outputImagePath);

            double mse = ImageUtils.calculateMSE(Imaging.getBufferedImage(new File("rhino.jpg")), outputImage);
            double psnr = ImageUtils.calculatePSNR(mse, 255);

            System.out.println("PSNR for HWT: " + psnr);
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ImageReadException e) {
            throw new RuntimeException(e);
        }
    }

    private static BufferedImage readImage(String filePath) throws IOException {
        return ImageIO.read(new File(filePath));
    }

    static void writeImage(BufferedImage image, String filePath) throws IOException {
        ImageIO.write(image, "jpg", new File(filePath)); // Change the format as needed
    }

    static double[][] convertTo2DUsingGetRGB(BufferedImage image) {
        int width = image.getWidth();
        int height = image.getHeight();
        double[][] result = new double[height][width];

        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                // In a grayscale image, red, green, and blue components are the same.
                // Thus, we can take any one of these values as the intensity.
                int intensity = new Color(image.getRGB(col, row)).getRed();
                result[row][col] = intensity;
            }
        }

        return result;
    }


    private static void applyHaarTransform(double[][] data, boolean isRowTransform) {
        int width = data[0].length;
        int height = data.length;

        // If isRowTransform is true, apply the transform to rows, else to columns
        for (int i = 0; i < (isRowTransform ? height : width); i++) {
            double[] temp;
            if (isRowTransform) {
                temp = data[i];
            } else {
                temp = new double[height];
                for (int j = 0; j < height; j++) {
                    temp[j] = data[j][i];
                }
            }

            // Apply Haar transform
            double[] transformed = haarTransform(temp);

            // Store the transformed data back
            if (isRowTransform) {
                data[i] = transformed;
            } else {
                for (int j = 0; j < height; j++) {
                    data[j][i] = transformed[j];
                }
            }
        }
    }

    private static double[] haarTransform(double[] arr) {
        int n = arr.length;
        double[] transformed = new double[n];

        while (n > 1) {
            n /= 2;
            for (int i = 0; i < n; i++) {
                transformed[i] = (arr[2 * i] + arr[2 * i + 1]) / 2;
                transformed[n + i] = (arr[2 * i] - arr[2 * i + 1]) / 2;
            }
            System.arraycopy(transformed, 0, arr, 0, n * 2);
        }

        return transformed;
    }

    private static void compressWaveletCoefficients(double[][] data, double threshold) {
        int width = data[0].length;
        int height = data.length;

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                // Thresholding: Set coefficients below a certain threshold to zero
                if (Math.abs(data[i][j]) < threshold) {
                    data[i][j] = 0;
                }
            }
        }
    }

    private static void applyInverseHaarTransform(double[][] data, boolean isRowTransform) {
        int width = data[0].length;
        int height = data.length;

        for (int i = 0; i < (isRowTransform ? height : width); i++) {
            double[] temp;
            if (isRowTransform) {
                temp = data[i];
            } else {
                temp = new double[height];
                for (int j = 0; j < height; j++) {
                    temp[j] = data[j][i];
                }
            }

            // Apply inverse Haar transform
            double[] transformed = inverseHaarTransform(temp);

            if (isRowTransform) {
                data[i] = transformed;
            } else {
                for (int j = 0; j < height; j++) {
                    data[j][i] = transformed[j];
                }
            }
        }
    }

    private static double[] inverseHaarTransform(double[] arr) {
        int length = arr.length;
        double[] transformed = Arrays.copyOf(arr, length);

        int n = 1;
        while (n < length) {
            double[] temp = Arrays.copyOf(transformed, length);
            for (int i = 0; i < n; i++) {
                int index1 = 2 * i;
                int index2 = 2 * i + 1;
                int index3 = n + i;

                if (index1 < length && index2 < length && index3 < length) {
                    transformed[index1] = temp[i] + temp[index3];
                    transformed[index2] = temp[i] - temp[index3];
                }
            }
            n *= 2;
        }

        return transformed;
    }


    private static BufferedImage reconstructImage(double[][] data) {
        int width = data[0].length;
        int height = data.length;
        BufferedImage reconstructedImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                // Clamp the value to the range 0-255
                int value = (int) Math.round(data[y][x]);
                value = Math.max(0, Math.min(value, 255)); // Ensure the value is between 0 and 255

                int rgb = new Color(value, value, value).getRGB();
                reconstructedImage.setRGB(x, y, rgb);
            }
        }
        return reconstructedImage;
    }


}

