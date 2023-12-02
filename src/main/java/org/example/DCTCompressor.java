package org.example;

import org.apache.commons.imaging.Imaging;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;

/*
Implementation of the Discrete Cosine Transform (DCT).

Step 1:
Convert image to block structure(as outlined in the task description, 8 by 8 blocks).

Step 2:
Apply DCT to each block.

Step 3:
Quantize the DCT coefficients.

Step 4:
Apply IDCT and reconstruct the image from quantized coefficients.
*/


public class DCTCompressor {

    // Constants
    private static final int BLOCK_SIZE = 8;
    private static final int QUANTIZATION_STEP = 150;

    public static void main(String[] args) {
        try {
            String inputImagePath = "rhino.jpg";
            String outputImagePath = "result_dct.jpg";

            // Step 1: Load image and convert to grayscale and block structure
            double[][] grayscaleImage = ImageUtils.loadImage(inputImagePath);

            // Step 2: Apply DCT to each block
            double[][][] dctTransformed = applyDCTToImage(grayscaleImage, BLOCK_SIZE);

            // Step 3: Quantize the DCT coefficients
            for (int i = 0; i < dctTransformed.length; i++) {
                for (int j = 0; j < dctTransformed[0].length; j++) {
                    // Flatten the block to 2D, quantize it, and then flatten back to 1D
                    double[][] twoDBlock = flattenTo2D(dctTransformed[i][j], BLOCK_SIZE);
                    double[][] quantizedBlock = quantizeBlock(twoDBlock, QUANTIZATION_STEP);
                    dctTransformed[i][j] = flattenToOneD(quantizedBlock, BLOCK_SIZE);
                }
            }

            // Step 4: Apply IDCT and reconstruct the image from quantized coefficients
            BufferedImage reconstructedImage = reconstructImageFromQuantized(dctTransformed, BLOCK_SIZE, QUANTIZATION_STEP);

            // Step 5: Save or display the reconstructed image
            ImageUtils.writeImage(reconstructedImage, outputImagePath);

            double mse = ImageUtils.calculateMSE(Imaging.getBufferedImage(new File("rhino.jpg")), reconstructedImage);
            double psnr = ImageUtils.calculatePSNR(mse, 255);

            System.out.println("PSNR for DCT: " + psnr);
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("An error occurred during the DCT compression process.");
        }
    }

    public static double[][][] applyDCTToImage(double[][] image, int blockSize) {
        int width = image.length;
        int height = image[0].length;
        double[][][] dctTransformed = new double[width / blockSize][height / blockSize][blockSize * blockSize];

        for (int i = 0; i < width; i += blockSize) {
            for (int j = 0; j < height; j += blockSize) {
                double[][] block = new double[blockSize][blockSize];
                for (int x = 0; x < blockSize; x++) {
                    System.arraycopy(image[i + x], j, block[x], 0, blockSize);
                }
                double[][] transformedBlock = applyDCT(block);
                // Flatten the 2D block into a 1D array for storage
                for (int x = 0; x < blockSize; x++) {
                    System.arraycopy(transformedBlock[x], 0, dctTransformed[i / blockSize][j / blockSize], x * blockSize, blockSize);
                }
            }
        }

        return dctTransformed;
    }


    public static double[][] applyDCT(double[][] block) {
        int size = block.length;
        double[][] dct = new double[size][size];

        for (int u = 0; u < size; u++) {
            for (int v = 0; v < size; v++) {
                double sum = 0.0;
                for (int i = 0; i < size; i++) {
                    for (int j = 0; j < size; j++) {
                        sum += block[i][j] * Math.cos(((2 * i + 1) * u * Math.PI) / (2 * size))
                                * Math.cos(((2 * j + 1) * v * Math.PI) / (2 * size));
                    }
                }
                sum *= ((u == 0) ? 1 / Math.sqrt(2) : 1) * ((v == 0) ? 1 / Math.sqrt(2) : 1);
                dct[u][v] = 0.25 * sum;
            }
        }

        return dct;
    }

    public static double[][] quantizeBlock(double[][] dctBlock, int quantizationStep) {
        int size = dctBlock.length;
        double[][] quantizedBlock = new double[size][size];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                quantizedBlock[i][j] = Math.round(dctBlock[i][j] / quantizationStep);
            }
        }

        return quantizedBlock;
    }

    public static double[][] applyIDCT(double[][] quantizedDctBlock) {
        int size = quantizedDctBlock.length;
        double[][] idctBlock = new double[size][size];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                double sum = 0.0;
                for (int u = 0; u < size; u++) {
                    for (int v = 0; v < size; v++) {
                        double C_u = u == 0 ? 1 / Math.sqrt(2) : 1;
                        double C_v = v == 0 ? 1 / Math.sqrt(2) : 1;
                        sum += C_u * C_v * quantizedDctBlock[u][v] *
                                Math.cos(((2 * i + 1) * u * Math.PI) / (2 * size)) *
                                Math.cos(((2 * j + 1) * v * Math.PI) / (2 * size));
                    }
                }
                idctBlock[i][j] = sum / 4;
            }
        }

        return idctBlock;
    }

    public static BufferedImage reconstructImageFromQuantized(double[][][] quantizedDctTransformed, int blockSize, int quantizationStep) {
        int width = quantizedDctTransformed.length * blockSize;
        int height = quantizedDctTransformed[0].length * blockSize;
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

        for (int i = 0; i < quantizedDctTransformed.length; i++) {
            for (int j = 0; j < quantizedDctTransformed[0].length; j++) {
                // De-quantization
                for (int u = 0; u < blockSize; u++) {
                    for (int v = 0; v < blockSize; v++) {
                        quantizedDctTransformed[i][j][u * blockSize + v] *= quantizationStep;
                    }
                }
                // Apply IDCT
                double[][] idctBlock = applyIDCT(flattenTo2D(quantizedDctTransformed[i][j], blockSize));
                // Reconstruct the image block by block
                for (int x = 0; x < blockSize; x++) {
                    for (int y = 0; y < blockSize; y++) {
                        int pixelValue = (int) (idctBlock[x][y] + 128);
                        pixelValue = Math.max(0, Math.min(255, pixelValue)); // Clamping to 0-255
                        Color color = new Color(pixelValue, pixelValue, pixelValue);
                        image.setRGB(i * blockSize + x, j * blockSize + y, color.getRGB());
                    }
                }
            }
        }
        return image;
    }

    public static double[][] flattenTo2D(double[] oneD, int blockSize) {
        double[][] twoD = new double[blockSize][blockSize];
        for (int i = 0; i < blockSize; i++) {
            System.arraycopy(oneD, i * blockSize, twoD[i], 0, blockSize);
        }
        return twoD;
    }

    public static double[] flattenToOneD(double[][] twoD, int blockSize) {
        double[] oneD = new double[blockSize * blockSize];
        for (int i = 0; i < blockSize; i++) {
            System.arraycopy(twoD[i], 0, oneD, i * blockSize, blockSize);
        }
        return oneD;
    }
}
