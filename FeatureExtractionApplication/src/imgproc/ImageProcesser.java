// Some par of this code: Descriptor creation, direction has influence from
// Code from http://aishack.in/tutorials/sift-scale-invariant-feature-transform-introduction/
// And https://github.com/aishack/sift
//Utkarsh Sinha Thank you for him.
// Also videos from youtube used
//https://www.youtube.com/watch?v=NPcMS49V5hg
//https://www.youtube.com/watch?v=oT9c_LlFBqs

package imgproc;

import static java.lang.Float.MAX_VALUE;
import static java.lang.Math.abs;
import static java.lang.Math.atan2;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.DataBufferByte;
import java.awt.image.WritableRaster;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

import javax.imageio.ImageIO;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.opencv.core.Core;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.MatOfKeyPoint;
import org.opencv.core.Size;
import org.opencv.features2d.Features2d;
import org.opencv.features2d.KeyPoint;
import org.opencv.highgui.Highgui;
import org.opencv.imgproc.Imgproc;

import m.LoggingConfig;

//
/**
 * Created by UDU on 17.07.2016.
 */
// https://medium.com/@aadimator/how-to-set-up-opencv-in-intellij-idea-6eb103c1d45c#.xotwy6mu0
//https://www.youtube.com/watch?v=123ZXYUJvu8
public class ImageProcesser {

	public static Logger logger;//
	public static int nrInstances = 0;
	static final int SIFT_FIXPT_SCALE = 1;
	public boolean firstPart = false;
	private int nrSmoothImgs = 5;
	private int nrOctaves = 4;
	private int nrIntervals = 2;
	private double M_PI = 3.1415926535897932384626433832795;
	// THESE ARE RECOMMENDED DEFAULT VARIABLES

	/// Parameters

	private final float img_scale = 1.f / (255 * SIFT_FIXPT_SCALE);
	private final float deriv_scale = img_scale * 0.5f;
	private final float second_deriv_scale = img_scale;
	private final float cross_deriv_scale = img_scale * 0.25f;

	private double Rbarrier = Double.MAX_VALUE;
	private int nrOfKeyPoints = 0;
	private int returnoctave = nrOctaves;
	private int returninterval = nrIntervals;
	private double sigma = 0.7;
	private double treshold = 2.5;
	private int BeginningBias = 0;// USED IN minMax, detectCorners
	private boolean MyMethod = false;
	private boolean GenerateComparator = true;
	private double defaultComparator = 0.03 * 255;
	private String imageName;
	private boolean shouldIprint = false;

	/// Tested these parameters as the BesT ,IMO.

	// Dastastructure
	ArrayList<double[]> subPixels = new <double[]>ArrayList();
	private double[][][] direction = new double[nrOctaves][nrIntervals][];
	private double[][][] magnitude = new double[nrOctaves][nrIntervals][];
	private double[][] sigmas = new double[nrOctaves][nrSmoothImgs];
	private ArrayList<ArrayList<BufferedImage>> dOgL = new <ArrayList<BufferedImage>>ArrayList();
	private ArrayList<ArrayList<BufferedImage>> kSmoothImgsbL = new <ArrayList<BufferedImage>>ArrayList();
	private ArrayList<ArrayList<BufferedImage>> keyPointsL = new <ArrayList<BufferedImage>>ArrayList();
	private ArrayList<ArrayList<ArrayList<MyKeyPoint>>> myKpL = new <ArrayList<ArrayList<MyKeyPoint>>>ArrayList();
	private ArrayList<ArrayList<ArrayList<KeyPoint>>> openKpL = new <ArrayList<ArrayList<KeyPoint>>>ArrayList();

	private BufferedImage originalImg;
	private BufferedImage forPrint;
	// private ArrayList<ArrayList<BufferedImage>> photoBuffer = new
	// ArrayList<ArrayList<BufferedImage>>();
	// Temporary list . goinna convert it to a proper thing later on; Yeah rigt
	// :D
	double[][][][] lSolutions = new double[nrOctaves][nrIntervals][][];

	public ImageProcesser(BufferedImage img, String pictureName) {
		// Logging Settings
		imageName = pictureName;
		if (nrInstances == 0) {
			logger = Logger.getLogger(getClass().getName());
			LoggingConfig.getlog(logger, 0);
		}
	}

	public void firstPart(BufferedImage img, String pictureName) {

		System.loadLibrary(Core.NATIVE_LIBRARY_NAME);
		originalImg = img;
		imageName = pictureName;
		myKpL.clear();

		for (int i = 0; i < nrOctaves; i++) {
			System.out.println(i);
			// Resize image
			BufferedImage gray = makeSmaller(img, i);
			// Resize image
			myKpL.add(new <ArrayList<MyKeyPoint>>ArrayList());
			kSmoothImgsbL.add(callGaussian1(gray, this.sigma * pow(2, i), i));// Also
																				// creates
																				// sigmas
			dOgL.add(substractImages(kSmoothImgsbL.get(i), i));
			keyPointsL.add(minMax(dOgL.get(i), gray.getWidth(), gray.getHeight(), i));// V繭tab
																						// dOgL-st
																						// i-nda
																						// oktavi
																						// pildid
			subPixel(i);
			for (int j = 0; j < nrIntervals; j++) {
				myKpL.get(i).add(detectCorners2(kSmoothImgsbL.get(i).get(j + this.BeginningBias),
						dOgL.get(i).get(j + this.BeginningBias), keyPointsL.get(i).get(j), this.treshold, i, j));
				// K�B L�B IGA OKTAAVI JA INTERVALLI
				// @parameter BufferedImage img - On this image creates a
				// specific R value and checks whether the KeyPoints are at all
				// good
				// I get this image from dOgL which containts 4 substracted
				// images for every octave. With bias 1 I do R procedures on
				// 1,2 . With bias 0,1 I dor R procedurse on 0,1
				// @parameter BufferedImage keyPoints - This is the image where
				// located minMax values are set to 255 and other values to 0
				// I get this image from keyPointsL which containt 2
				// (intervalNr) images for every octave.
			}
		}

		this.myKpL = geneRateHistogram(this.myKpL, this.direction, this.magnitude);
		this.openKpL = extraxtKeyPoints(this.myKpL);
		this.Rbarrier = reduceKP(this.openKpL);
		extractDescriptor(this.myKpL);

		// System.out.println(this.myKpL.get(0).get(1).get(0).getAngle_deg().get(0));
		int saveOctaves[] = { 0 };
		int saveScales[] = { 0, 1 };

		saveKeyPointImage(this.openKpL, false, saveOctaves, saveScales);
		System.out.println("Subpiksel size " + subPixels.size());
	}

	public double gtBtwnS(byte[] img, double pos1) {
		return 0.5 * (img[(int) (pos1 - 0.5)] + img[(int) (pos1 + 0.5)]);
	}

	public double gtBtwnV(byte[] img, double pos1, int width) {
		return 0.5 * img[(int) (1)];
	}

	public void extractDescriptor(ArrayList<ArrayList<ArrayList<MyKeyPoint>>> in) {
		// Interpolation
		double[][][] interPolatMag = new double[nrOctaves][][];
		double[][][] interPolatOrie = new double[nrOctaves][][];
		for (int i = 0; i < nrOctaves; i++) {
			int width = kSmoothImgsbL.get((int) i).get(0).getWidth();
			int height = kSmoothImgsbL.get((int) i).get(0).getHeight();
			int size = (width + 1) * (height + 1);
			interPolatMag[i] = new double[nrIntervals][size];
			interPolatOrie[i] = new double[nrIntervals][size];
			for (int j = 0; j < nrIntervals; j++) {

				BufferedImage temp = kSmoothImgsbL.get(i).get(j);
				byte[] data = ((DataBufferByte) temp.getRaster().getDataBuffer()).getData();

				for (double ii = 1.5; ii < height - 1.5; ii++) {
					for (double jj = 1.5; jj < width - 1.5; jj++) {
						double up15 = (ii - 1.5) * width + jj;
						double down15 = (ii + 1.5) * width + jj;
						double right15 = ii * width + jj + 1.5;
						double left15 = ii * width + jj - 1.5;
						double up05 = (ii - 0.5) * width + jj;
						double down05 = (ii + 0.5) * width + jj;
						double right05 = ii * width + jj + 0.5;
						double left05 = ii * width + jj - 0.5;
						double dx = data[(int) right15] + data[(int) right05] - data[(int) left15] - data[(int) left05];
						double dy = data[(int) up15] + data[(int) up05] - data[(int) down15] - data[(int) down05];
						int iii = (int) ii + 1;
						int jjj = (int) jj + 1;
						interPolatMag[i][j][iii * width + jjj] = sqrt((dx * dx + dy * dy));
						interPolatOrie[i][j][iii * width + jjj] = atan2(dy, dx);
					}
				}
				BufferedImage forSaving1 = new BufferedImage(width + 1, height + 1, BufferedImage.TYPE_BYTE_GRAY);
				BufferedImage forSaving2 = new BufferedImage(width + 1, height + 1, BufferedImage.TYPE_BYTE_GRAY);
				byte[] forSave1 = ((DataBufferByte) forSaving1.getRaster().getDataBuffer()).getData();
				byte[] forSave2 = ((DataBufferByte) forSaving2.getRaster().getDataBuffer()).getData();

				double maxMag = interPolatMag[i][j][0];
				double maxOrie = interPolatOrie[i][j][0];
				double minMag = interPolatMag[i][j][0];
				double minOrie = interPolatOrie[i][j][0];

				for (int l = 0; l < size; l++) {

					maxMag = maxMag > interPolatMag[i][j][l] ? maxMag : interPolatMag[i][j][l];
					maxOrie = maxOrie > interPolatOrie[i][j][l] ? maxOrie : interPolatOrie[i][j][l];
					minMag = minMag < interPolatMag[i][j][l] ? minMag : interPolatMag[i][j][l];
					minOrie = minOrie > interPolatOrie[i][j][l] ? minOrie : interPolatOrie[i][j][l];
				}

				for (int l = 0; l < size; l++) {
					forSave1[l] = (byte) (interPolatMag[i][j][l]);
					forSave2[l] = (byte) (interPolatOrie[i][j][l]);
				}

				saveImage(forSaving1, "MAG_O" + i + "I" + j, "InterpolTest/");
				saveImage(forSaving2, "ORIE_O" + i + "I" + j, "InterpolTest/");

			}
		}
	}

	public ArrayList<MyKeyPoint> detectCorners2(BufferedImage forGradient, BufferedImage img, BufferedImage keyPoints,
			double threshold, int octave, int scale) {

		BufferedImage img_copy1 = deepCopy(img);
		BufferedImage img_copy2 = deepCopy(img);

		ArrayList<MyKeyPoint> in = new <MyKeyPoint>ArrayList();
		int width = img.getWidth();
		int height = img.getHeight();
		byte[] f = ((DataBufferByte) img_copy1.getRaster().getDataBuffer()).getData();
		byte[] keypoint = ((DataBufferByte) keyPoints.getRaster().getDataBuffer()).getData();
		double R[] = getR(f, keypoint, img.getHeight(), img.getWidth(), threshold); // Sets
																					// some
																					// keypoint
																					// values
																					// to
																					// 0
		double temp[][] = gradientImage(forGradient, img.getWidth(), img.getHeight());

		// Adds an

		this.magnitude[octave][scale] = temp[0];
		this.direction[octave][scale] = temp[1];

		for (int y = 1; y < height - 1; y++) {
			for (int x = 1; x < width - 1; x++) {
				// logger.info("xy" + x +" "+y);
				double val = R[y * width + x];
				if (keypoint[y * width + x] != 0 && checkVal(R, y, x, val, width)) {
					MyKeyPoint kp = new MyKeyPoint(x, y, octave, scale, R[y * width + x]);
					// kp.setDegree((float)
					// direction[octave][scale][y*width+x]);
					kp.setSigma((float) sigmas[octave][scale]);
					in.add(kp);
				} else {
					keypoint[y * img.getWidth() + x] = 0;
				}
			}
		}

		saveImage(keyPoints, "Keypoints_o" + octave + "s" + scale, "Corner/");
		saveImage(img, "Input_o" + octave + "s" + scale, "Corner/");
		return in;
	}

	public ArrayList<ArrayList<ArrayList<KeyPoint>>> extraxtKeyPoints(ArrayList<ArrayList<ArrayList<MyKeyPoint>>> m) {
		ArrayList<KeyPoint> temp = new <KeyPoint>ArrayList();
		ArrayList<ArrayList<ArrayList<KeyPoint>>> returnValue = new <ArrayList<ArrayList<KeyPoint>>>ArrayList();
		int rejected = 0;
		for (int i = 0; i < nrOctaves; i++) {
			returnValue.add(new <ArrayList<KeyPoint>>ArrayList());
			for (int j = 0; j < nrIntervals; j++) {
				returnValue.get(i).add(new <KeyPoint>ArrayList());
				for (MyKeyPoint kp : m.get(i).get(j)) {
					if (kp.getListsize() != 0) {
						temp = kp.getOpenCvKp();
						// returnValue.get(i).get(j).add(kp.returnSimple());
						for (KeyPoint l : temp) {
							nrOfKeyPoints++;
							returnValue.get(i).get(j).add(l);
						}
					} else {
						rejected++;
					}

				}
			}
		}
		logger.fine(this.getClass().getName() + "Rejected" + rejected);
		return returnValue;
	}

	public ArrayList<ArrayList<ArrayList<MyKeyPoint>>> geneRateHistogram(ArrayList<ArrayList<ArrayList<MyKeyPoint>>> in,
			double[][][] drection, double[][][] magnitude) {
		int failCount = 0;
		int goodCount = 0;
		double[] histogram = new double[36];
		for (int i = 0; i < nrOctaves; i++) {
			int scale = (int) pow(2.0, (double) i);
			for (int j = 0; j < nrIntervals; j++) {
				double abs_sigma = sigmas[i][j];
				int width = this.keyPointsL.get(i).get(j).getWidth();
				int height = this.keyPointsL.get(i).get(j).getHeight();
				gaussianSmoothing((magnitude[i][j]), 1.5 * abs_sigma, this.dOgL.get(i).get(j).getWidth(),
						this.dOgL.get(i).get(j).getHeight());
				int hfsz = createGaussianMask3(1.5 * abs_sigma).length / 2;
				// double[] mask = createGaussianMaskM(1.5*abs_sigma);
				int size = in.get(i).get(j).size();

				if (size != 0) {
					for (int l = 0; l < size; l++) {
						MyKeyPoint kp = null;
						try {
							kp = in.get(i).get(j).get(l);
						} catch (Exception kk) {
							kk.printStackTrace();
							logger.severe("Cant load keypoint for histogram");
						}

						for (int k = 0; k < 36; k++) {
							histogram[k] = 0.0;
						}
						if (kp != null) {
							kp.addSize(hfsz);
							int x = kp.getX();
							int y = kp.getY();
							// Check the surrounding to create the histogram;
							for (int kk = -hfsz; kk <= hfsz; kk++) {
								for (int tt = -hfsz; tt <= hfsz; tt++) {
									if ((x + tt < 0 || x + tt >= width || y + kk < 0 || y + kk >= height) == false) {

										double sampleOrient = drection[i][j][(y + kk) * width + x + tt];
										if (sampleOrient <= -M_PI || sampleOrient >= M_PI)
											logger.severe("Bad Corner" + sampleOrient);

										sampleOrient += M_PI;
										double sampleOrientDegrees = abs(sampleOrient * 180 / M_PI);
										// hist_orient[(int)sampleOrientDegrees
										// / (360/NUM_BINS)] +=
										// cvGetReal2D(imgWeight, yi+tt, xi+kk);
										// System.out.println(magnitude[i][j][(y
										// + kk) * width + x + tt]);
										histogram[(int) (sampleOrientDegrees
												/ (10))] += (magnitude[i][j][(y + kk) * width + x + tt]) / 256;
									}
								}
							}

							double max_peak = histogram[0];
							int max_peak_index = 0;
							for (int k = 0; k < 36; k++) {
								if (histogram[k] > max_peak) {
									max_peak = histogram[k];
									max_peak_index = k;
								}
							}

							// System.out.println(max_peak);

							// Next we fit a downward parabola aound
							// these three points for better accuracy

							// A downward parabola has the general form
							//
							// y = a * x^2 + bx + c
							// Now the three equations stem from the three
							// points
							// (x1,y1) (x2,y2) (x3.y3) are
							//
							// y1 = a * x1^2 + b * x1 + c
							// y2 = a * x2^2 + b * x2 + c
							// y3 = a * x3^2 + b * x3 + c
							//
							// in Matrix notation, this is y = Xb, where
							// y = (y1 y2 y3)' b = (a b c)' and
							//
							// x1^2 x1 1
							// X = x2^2 x2 1
							// x3^2 x3 1
							//
							// OK, we need to solve this equation for b
							// this is done by inverse the matrix X
							//
							// b = inv(X) Y
							// Here one or many keypoint corners is added

							if (max_peak != 0) {
								ArrayList<Double> orien;
								ArrayList<Double> mag;
								for (int k = 0; k < 36; k++) {
									// IF we have a good peak, then we add a
									// direction to our Keypoint
									if (histogram[k] > 0.8 * max_peak) {

										// Three points. (x2,y2) is the peak and
										// (x1,y1)
										// and (x3,y3) are the neigbours to the
										// left and right.
										// If the peak occurs at the extreme
										// left, the "left
										// neighbour" is equal to the right
										// most. Similarly for
										// the other case (peak is rightmost)
										double x1 = k - 1;
										double y1;
										double x2 = k;
										double y2 = histogram[k];
										double x3 = k + 1;
										double y3;

										if (k == 0) {
											y1 = histogram[36 - 1];
											y3 = histogram[1];
										} else if (k == 35) {
											y1 = histogram[36 - 1];
											y3 = histogram[0];
										} else {
											y1 = histogram[k - 1];
											y3 = histogram[k + 1];
										}

										double[][] matrixData = { { x1 * x1, x1, 1 }, { x2 * x2, x2, 1 },
												{ x3 * x3, x3, 1 } };
										RealMatrix m = MatrixUtils.createRealMatrix(matrixData);
										double[][] y_matrix = { { y1 }, { y2 }, { y3 } };
										RealMatrix y_real = MatrixUtils.createRealMatrix(y_matrix);
										RealMatrix mInverse = new LUDecomposition(m).getSolver().getInverse();
										RealMatrix result = mInverse.multiply(y_real);
										double a = result.getEntry(0, 0);
										double b = result.getEntry(1, 0);
										double c = result.getEntry(2, 0);

										double parabola_top = -b / 2 * a;
										if (abs(parabola_top) > 2 * 36) {
											failCount++;
											// logger.severe("Parabola top is a
											// very wrong number"+histogram[k]);
											parabola_top = x2;
										} else if (parabola_top < 0) {
											failCount++;
											while (parabola_top < 0) {
												parabola_top += 36;
											}
										} else if (parabola_top > 36) {
											failCount++;
											while (parabola_top > 36) {
												parabola_top -= 36;
											}
										} else {
											goodCount++;
										}
										// Normalize
										double parabola_topN = M_PI * parabola_top / 36;// Scaling
																						// min
																						// 0
																						// max
																						// 2*M_PI
										if (!(parabola_topN >= 0 && parabola_topN < 2 * M_PI)) {
											logger.severe("Normalization of corner went wrong" + parabola_top);
										}
										/*
										 * parabola_topN -= M_PI;
										 * if(!(parabola_topN>=-M_PI &&
										 * parabola_topN<M_PI)){ logger.
										 * severe("Normalization of corner went wrong"
										 * +parabola_top); }
										 */

										kp.addAngle((float) (parabola_topN / (2 * M_PI) * 360));
										kp.addAngleMag(histogram[k]);

									}
								}
							}
						}
					}
				}
			}
		}
		logger.fine("GoodCount" + goodCount);
		logger.fine("failCount" + failCount);
		return in;
	}

	public boolean checkVal(double[] R, int y, int x, double val, int width) {
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				if (val > R[(y + i) * width + x + j]) {
					return false;
				}
			}
		}
		return true;
	}

	/// R is very well explained in the lecture slides. Basically. What it does
	/// is approximating eigenvalue relative size to each other. We are looking
	/// for the casees wher both eigenvalues are large.

	public double[] getR(byte[] f, byte keypoint[], int height, int width, double threshold) {

		double trace2 = 0;
		double determinant = 0;
		int allkeypoints = 0;
		int notAccepted = 0;
		double[] R = new double[f.length];
		for (int x = 1; x < height - 1; x++) {
			for (int y = 1; y < width - 1; y++) {
				if (keypoint[x * width + y] != 0) {
					allkeypoints++;
					double fxx = (((int) f[(x + 1) * width + y] & 0xFF) + ((int) f[(x - 1) * width + y] & 0xFF)
							- 2 * f[x * width + y]) * second_deriv_scale;
					double fyy = (((int) f[(x * width) + y + 1] & 0xFF) + ((int) f[(x * width) + y - 1] & 0xFF)
							- 2 * f[x * width + y]) * second_deriv_scale;
					double fx_y = (((int) f[(x + 1) * width + y + 1] & 0xFF) + ((int) f[(x - 1) * width + y - 1] & 0xFF)
							- ((int) f[(x - 1) * width + y + 1] & 0xFF) - ((int) f[(x - 1) * width + y + 1] & 0xFF))
							* cross_deriv_scale;

					trace2 = (fxx + fyy) * (fxx + fyy);// (fx2[x*width+y]+fy2[x*width+y])*(fx2[x*width+y]+fy2[x*width+y]);
					determinant = fxx * fyy - fx_y * fx_y;// fx2[x*width+y]*fy2[x*width+y]
															// -
															// fxy[x*width+y]*fxy[x*width+y];
					R[x * width + y] = abs((trace2) / (determinant + 0.000000000001));
					if (determinant > 0 && abs(R[x * width + y]) < (threshold + 1) * (threshold + 1) / (threshold)) {
					} else {
						keypoint[x * width + y] = 0;
						R[x * width + y] = threshold * threshold;
						notAccepted++;
					}
				} else {
					R[x * width + y] = threshold * threshold;
				}
			}
		}

		logger.fine("All keyPoints:" + allkeypoints);
		logger.fine("Dismissed keyPoints by R:" + notAccepted);

		return R;
	}

	// EXTREMA FUNCTIONS

	public ArrayList<BufferedImage> minMax(ArrayList<BufferedImage> list, int width, int height, int octave) {

		// IM GONNA A CREATE A NEW IMAGE LIST FOR RETURNING KEYPOINTS AND A LIST
		// IM GONNA RETURN THESE

		ArrayList<BufferedImage> returnValue = new <BufferedImage>ArrayList();
		BufferedImage temp1 = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
		BufferedImage temp2 = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
		returnValue.add(temp1);
		returnValue.add(temp2);

		// HAVE TO CONVERT THEM TO BYTE TO ADD THE COMPARISON RESULTS TO OUTPUT

		byte[] result1 = ((DataBufferByte) temp1.getRaster().getDataBuffer()).getData();
		byte[] result2 = ((DataBufferByte) temp2.getRaster().getDataBuffer()).getData();
		byte[][] minMax = new byte[2][];
		minMax[0] = result1;
		minMax[1] = result2;

		// INITIALIZE A NEW BYTE AND DOUBLE LIST BASED ON INPUT FOR COMPARING
		// Calculate average to form a threshold

		double totalSum = 0.0;
		int[] arrai = new int[256];
		byte[][] dOgb = new byte[list.size()][];
		double[][] dOgd = new double[list.size()][];
		for (int i = 0; i < list.size(); i++) {
			dOgb[i] = ((DataBufferByte) list.get(i).getRaster().getDataBuffer()).getData();
			dOgd[i] = new double[dOgb[i].length];
			for (int j = 0; j < dOgb[i].length; j++) {
				dOgd[i][j] = (double) (dOgb[i][j] & 0xFF);
				arrai[dOgb[i][j] & 0xFF]++;
				totalSum++;
			}
		}

		double comparator = this.defaultComparator;
		if (this.GenerateComparator) {
			for (int i = 1; i < 256; i++) {
				arrai[i] += arrai[i - 1];
				if (arrai[i] / totalSum > 0.90) {
					comparator = (double) i;
					break;
				}
			}
		}

		// logger.info("Comparator " + comparator);

		// SITUATION: I HAVE A LIST BASED ON dOg. AKA 4 IMAGES BASED ON
		// SUBSTRACTION. MY PURPOSE IS TO CREATE 2 IMAGES BASED
		// ON 2 3 ELEMENT SUBSETS OF THESE 4 IMAGES ({0,1,2},{1,2,3}). CHOOSING
		// MAXIMA AND MINIMA
		double max = 0;
		for (int i = 0; i < 2; i++) {// i = 0,1
			for (int y = 1; y < height - 1; y++) {
				for (int x = 1; x < width - 1; x++) {
					if (this.MyMethod) {
						max = getMax(dOgd[i][y * width + x], dOgb[i + 1][y * width + x], dOgb[i + 2][y * width + x]);// CHOOSING
																														// THE
																														// SUBSETSmax
																														// =
																														// dOgb[i+1][y
																														// *
																														// width
																														// +
																														// x];
					} else {
						max = dOgb[i + this.BeginningBias][y * width + x];
					}
					// MAX BASED ON 3 IMG IN COORD X AND Y
					if (max > comparator) {
						if ((check(dOgd, y, x, max, width, i + 1)) == true) { // MIDDLEPOSITION
																				// EITHER
																				// 1
																				// OR
																				// 2
							minMax[i][y * width + x] = (byte) 255;
						} else if ((checkf(dOgd, y, x, max, width, i + 1)) == true) {
							minMax[i][y * width + x] = (byte) 255;
						} else {
							minMax[i][y * width + x] = 0;
						}
					}
				}
			}
		}

		String m1 = "m0+" + octave + "Size" + height + "X" + width;
		String m2 = "m1" + octave + "Size" + height + "X" + width;
		saveImage(temp1, m1, "MinMaxTest/");
		saveImage(temp2, m2, "MinMaxTest/");
		return returnValue;

	}

	void subPixel(int octave) {

		double[] ad;
		int i = octave;
		ArrayList<BufferedImage> temp = dOgL.get(i);
		// System.out.println("DIS IS i"+i);
		int height = temp.get(0).getHeight();
		int width = temp.get(0).getWidth();
		for (int s = 0; s < nrIntervals; s++) {
			byte[] upaluu = ((DataBufferByte) keyPointsL.get(i).get(s).getRaster().getDataBuffer()).getData();
			// System.out.println("s" +s);
			for (int y = 0; y < height - 5; y++) {
				for (int x = 0; x < width - 5; x++) {
					if (upaluu[y * width + x] != 0) {
						int lay_temp = s;
						double y_temp = y;
						double x_temp = x;
						for (int m = 0; m < 5; m++) {
							byte[] down = ((DataBufferByte) temp.get(lay_temp).getRaster().getDataBuffer()).getData();
							byte[] middle = ((DataBufferByte) temp.get(lay_temp + 1).getRaster().getDataBuffer())
									.getData();
							byte[] up = ((DataBufferByte) temp.get(lay_temp + 2).getRaster().getDataBuffer()).getData();
							// TULEB SAADA K�TE x+1,x-1,y+1,y-1,z+1,z-1
							// float v2 = (float)img.at<sift_wt>(r, c)*2;
							int iwj = (int) (y_temp * width + x_temp);
							int vertical = width;
							double M_C = middle[iwj];
							double M_R = middle[iwj + 1];
							double M_L = middle[iwj - 1];
							double M_A = middle[iwj + vertical];
							double M_U = middle[iwj - vertical];
							double M_DSum1 = middle[iwj + vertical + 1] + middle[iwj - vertical - 1];
							double M_Dsum2 = middle[iwj - vertical + 1] + middle[iwj + vertical - 1];
							double dxx = (M_R + M_L - 2 * M_C) * second_deriv_scale;
							double dyy = (M_A + M_U - 2 * M_C) * second_deriv_scale;
							double dxy = (M_DSum1 - M_Dsum2) * cross_deriv_scale;
							double dss = (up[iwj] + down[iwj] - 2 * middle[iwj]) * second_deriv_scale;
							double dxs = (up[iwj + 1] - up[iwj - 1] - down[iwj + 1] + down[iwj - 1])
									* cross_deriv_scale;
							double dys = (up[iwj + vertical] - up[iwj - vertical] - down[iwj + vertical]
									+ down[iwj - vertical]) * cross_deriv_scale;
							double[][] matrixData = { { dxx, dxy, dxs }, { dxy, dyy, dys }, { dxs, dys, dss } };
							double[][] dp = { { (M_R - M_L) }, { (M_A - M_U) }, { (up[iwj] - down[iwj]) } };
							RealMatrix md = MatrixUtils.createRealMatrix(matrixData);
							RealMatrix dpd = MatrixUtils.createRealMatrix(dp);
							// System.out.println(md);
							// System.out.println(dpd);
							try {
								RealMatrix mdInverse = new LUDecomposition(md).getSolver().getInverse();
								RealMatrix result = mdInverse.multiply(dpd);
								double xi = -result.getEntry(2, 0);
								double xr = -result.getEntry(1, 0);
								double xc = -result.getEntry(2, 0);
								// System.out.println(xi+" "+xr+" "+xc);
								ad = new double[] { i, s, y, x, xi, xr, xc };
								if (abs(xi) < 0.5f && abs(xr) < 0.5f && abs(xc) < 0.5f) {
									subPixels.add(ad);
									break;
								}
								if (abs(xi) > (float) (MAX_VALUE / 3) || abs(xr) > (float) (MAX_VALUE / 3)
										|| abs(xc) > (float) (MAX_VALUE / 3)) {
									break;
								} // DONOTHING AND FO GOR NEXT X,Y
								x_temp += (int) xc;
								y_temp += (int) xr;
								lay_temp += (int) xi;
								if (lay_temp < 0 || lay_temp >= nrIntervals || x_temp < 5 || x_temp >= width - 5
										|| y_temp < 5 || y_temp >= height - 5) {
									break;
								} // DONOTHING AND FO GOR NEXT X,Y
							} catch (SingularMatrixException k) {
								// k.printStackTrace();
								// System.out.println("Singular matrix");
								break;
							}

						}

					}
					// N� TULEB AND KEYPOINDILE UUS V�TUS
				}
			}
		}

	}

	// I WOULD SAY THIS IS THE MOST IMPORTANT FUNCTION FOR MINMAX
	public boolean check(double[][] dOgd, int y, int x, double max, int width, int middleposition) {
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				for (int l = -1; l <= 1; l++) {
					if (max < dOgd[middleposition + i][(y + j) * width + x + l]) {
						return false;
					}
				}
			}
		}
		return true;
	}

	public boolean checkf(double[][] dOgd, int y, int x, double max, int width, int middleposition) {
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				for (int l = -1; l <= 1; l++) {
					if (max > dOgd[middleposition + i][(y + j) * width + x + l]) {
						return false;
					}
				}
			}
		}
		return true;
	}

	public double[][] gradientImage(BufferedImage image, int width, int height) {

		double[] img1 = new double[width * height];
		double[] img2 = new double[width * height];
		double[][] returnVal = { img1, img2 };

		byte[] original = ((DataBufferByte) image.getRaster().getDataBuffer()).getData();
		int[] filter = new int[] { -1, 0, 1 };
		int a = filter.length / 2;

		for (int x = a; x < height - a; x++) {
			for (int y = a; y < width - a; y++) {
				double sum_x = 0.0;
				double sum_y = 0.0;
				sum_x += (original[(x + 1) * width + y] & 0xFF) - (original[(x - 1) * width + y] & 0xFF);
				sum_y += (original[(x) * width + y + 1] & 0xFF) - (original[(x) * width + y - 1] & 0xFF);
				returnVal[0][x * width + y] = (Math.sqrt((sum_x * sum_x + sum_y * sum_y)));
				returnVal[1][x * width + y] = Math.atan(sum_x / (sum_y + 0.0000000001));
			}
		}

		// IMAGE TEST SAVE
		BufferedImage mag = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
		BufferedImage dir = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
		byte[] magb = ((DataBufferByte) mag.getRaster().getDataBuffer()).getData();
		byte[] dirb = ((DataBufferByte) dir.getRaster().getDataBuffer()).getData();

		for (int i = 0; i < width * height; i++) {
			dirb[i] = (byte) (returnVal[1][i] / M_PI * 256);
			magb[i] = (byte) (30 * returnVal[0][i]);
		}
		saveImage(mag, "Magnitude" + height + "X" + width, "GradMagnTest/");
		saveImage(dir, "Direction" + height + "X" + width, "GradMagnTest/");
		//

		return returnVal;
	}

	public double getMax(double a, double b, double c) {
		double maxI = abs(a) > abs(b) ? a : b;
		return abs(maxI) > abs(c) ? maxI : c;
	}

	// SUBSTRACT FUNCTIONS

	private ArrayList<BufferedImage> substractImages(ArrayList<BufferedImage> kSmoothI, int octave) {

		ArrayList<BufferedImage> temp = new <BufferedImage>ArrayList();
		for (int i = 1; i < nrSmoothImgs; i++) {
			temp.add(substract(kSmoothI.get(i), kSmoothI.get(i - 1)));
		}

		for (int i = 0; i < nrSmoothImgs - 1; i++) {
			BufferedImage m = temp.get(i);
			saveImage(m, "SubsNr_" + i + " " + octave + "" + m.getHeight() + "X" + m.getWidth(), "DogTest/");
		}
		return temp;
	}

	private BufferedImage substract(BufferedImage img1, BufferedImage img2) {
		int size = img1.getHeight() * img1.getWidth();
		BufferedImage temp = new BufferedImage(img1.getWidth(), img1.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
		byte[] result = ((DataBufferByte) temp.getRaster().getDataBuffer()).getData();
		byte[] imgb1 = ((DataBufferByte) img1.getRaster().getDataBuffer()).getData();
		byte[] imgb2 = ((DataBufferByte) img2.getRaster().getDataBuffer()).getData();
		for (int i = 0; i < size & i < size; i++) {
			result[i] = (byte) abs(((imgb1[i] & 0xFF) - (imgb2[i] & 0xFF)));
		}
		return temp;
	}

	// http://www.javaworld.com/article/2077578/learn-java/java-tip-76--an-alternative-to-the-deep-copy-technique.html
	//
	public BufferedImage deepCopy(BufferedImage bi) {
		ColorModel cm = bi.getColorModel();
		boolean isAlphaPremultiplied = cm.isAlphaPremultiplied();
		WritableRaster raster = bi.copyData(bi.getRaster().createCompatibleWritableRaster());
		return new BufferedImage(cm, raster, isAlphaPremultiplied, null);
	}

	// SMOOTHER SCALING FUNCTIONS

	private void gaussianSmoothing(double[] f, double sigma, int width, int height) {
		double[] gMask = createGaussianMask3(sigma);
		double[] img = new double[f.length];
		convolve4(f, width, height, gMask, gMask.length / 2, 0);
		convolve4(f, width, height, gMask, 0, gMask.length / 2);
		// convert back to image format

	}

	private void convolve4(double[] f, int width, int height, double[] mask, int a, int b) {

		double[] f2 = new double[f.length];
		System.arraycopy(f, 0, f2, 0, f.length);
		double sum = 0.0;
		int filterSize = Math.min(2 * a + 1, 2 * b + 1);
		for (int x = a; x < height - a; x++) {
			for (int y = b; y < width - b; y++) {
				sum = 0.0;
				for (int s = -a; s <= +a; s++) {
					for (int t = -b; t <= +b; t++) {
						sum += mask[(s + a) * filterSize + t + b] * (f2[(x + s) * width + y + t]);
					}
					f[x * width + y] = sum;
				}
			}
		}
	}

	public ArrayList<BufferedImage> callGaussian1(BufferedImage img, double sigma, int octave) {

		logger.entering(getClass().getName(), "callSmoothing");
		int width = img.getWidth();
		int height = img.getHeight();
		ArrayList<BufferedImage> temp = new <BufferedImage>ArrayList();
		if (img != null) {
			BufferedImage imgCopy;// = new BufferedImage(width, height,
									// BufferedImage.TYPE_BYTE_GRAY);
			imgCopy = deepCopy(img);
			byte[] imgData = ((DataBufferByte) imgCopy.getRaster().getDataBuffer()).getData(); // new
																								// byte[imgDataOrg.length];
			if (imgData != null) {
				for (int i = 0; i < nrSmoothImgs; i++) {
					logger.fine((getClass().getName() + "For loop nr:" + i).toString());
					if (imgData != null) {
						// System.out.printnln()
						sigmas[octave][i] = (sigma) * pow(sqrt(2), i + 1);
						logger.fine("Sigma: " + (sigma) * pow(sqrt(2), i + 1));
						gaussianSmoothing2(imgData, (sigma) * pow(sqrt(2), (i + 1)), img.getWidth(), img.getHeight());
						temp.add(deepCopy(imgCopy));
					}
				}
			} else {
				logger.severe("Did not load image data");
			}
		}

		int cnt = 0;
		logger.exiting(getClass().getName(), "callSmoothing");
		for (BufferedImage i : temp) {
			saveImage(i, "oct" + octave + " scal" + cnt + "" + img.getWidth(), "SmoothnessTest/");
			cnt++;
		}
		return temp;
	}

	// Intellijj offered this to me automatically. This is a good idea for
	// default parameters.
	/*
	 * private void gaussianSmoothing(int width) { gaussianSmoothing(,, width,);
	 * }
	 */

	private void gaussianSmoothing2(byte[] f, double sigma, int width, int height) {
		double[] gMask = createGaussianMask3(sigma);
		double[] img = new double[f.length];
		// Convert to double for better playing with numbers
		for (int i = 0; i < f.length; i++) {
			img[i] = f[i] & 0xFF;
		}
		convolve4(img, width, height, gMask, gMask.length / 2, 0);
		convolve4(img, width, height, gMask, 0, gMask.length / 2);
		// convert back to image format
		for (int i = 0; i < f.length; i++) {
			f[i] = (byte) img[i];
		}

	}

	double[] getDerivative(int x, int y, BufferedImage img) {

		if (x > 2 && y > 2 && x < (img.getWidth() - 2) && y < (img.getHeight() - 2)) {
			byte[] f = ((DataBufferByte) img.getRaster().getDataBuffer()).getData(); // new
																						// byte[imgDataOrg.length];
			int width = img.getWidth();
			double fxx = ((int) f[(x + 1) * width + y] & 0xFF) + ((int) f[(x - 1) * width + y] & 0xFF)
					- 2 * f[x * width + y];
			double fyy = ((int) f[(x * width) + y + 1] & 0xFF) - ((int) f[(x * width) + y - 1] & 0xFF)
					- 2 * f[x * width + y];
			double fx_y = ((int) f[(x - 1) * width + y + 1] & 0xFF) + ((int) f[(x - 1) * width + y + 1] & 0xFF)
					- ((int) f[(x + 1) * width + y + 1] & 0xFF) - ((int) f[(x - 1) * width + y - 1] & 0xFF);
			double trace2 = (fxx + fyy) * (fxx + fyy);// (fx2[x*width+y]+fy2[x*width+y])*(fx2[x*width+y]+fy2[x*width+y]);
			double determinant = fxx * fyy - fx_y * fx_y;// fx2[x*width+y]*fy2[x*width+y]
															// -
															// fxy[x*width+y]*fxy[x*width+y];
			double R = abs((trace2) / (determinant + 0.000000000001));
			return new double[] { fxx, fyy, fx_y, trace2, determinant, R };
		}
		return new double[] { 0, 0, 0, 0, 0, 0 };

	}

	public double[] createGaussianMask3(double sigma) {
		ArrayList<Double> gaussMask = new ArrayList();
		double gaussianDieOff = 0.01;
		double d;
		for (double i = 0; true; i++) {
			d = Math.exp(-i * i / (2.0 * sigma * sigma));
			if (d < gaussianDieOff)
				break;
			gaussMask.add(d);
		}

		int gSize = gaussMask.size();
		int l = ((gaussMask.size() - 1) * 2) + 1;
		double[] fullMask = new double[l];
		double sumOfGaussian = 0.0;
		for (int i = 0; i < l; i++) {
			int k = gSize - (i + 1);
			double b;
			b = (double) k + 0.1;
			fullMask[i] = Math.abs(gaussMask.get(Math.abs(k)));
			sumOfGaussian += fullMask[i];
		}

		logger.finest((getClass().getName() + "Gaussian mask").toString());
		String a = "";
		for (int i = 0; i < fullMask.length; i++) {
			fullMask[i] = fullMask[i] / sumOfGaussian;
			a += (" " + fullMask[i]).toString();

		}
		logger.finest(a);
		return fullMask;
	}

	public void saveImage(BufferedImage img, String addition, String folder) {
		String name = "Result" + addition;
		boolean notOk = true;
		File output = new File("./Tests/" + folder + name + ".png");
		try {
			ImageIO.write(img, "png", output);
		} catch (Exception ell) {
			System.out.println("Couldnt save image");
		}

	}

	public double reduceKP(ArrayList<ArrayList<ArrayList<KeyPoint>>> in) {
		// MatOfKeyPoint keypoints_s = new MatOfKeyPoint();
		double[] Rlist = new double[nrOfKeyPoints];
		int cnt = 0;
		if (true) {
			for (int i = 0; i < nrOctaves; i++) {// OCTAVE LEVEL
				for (int j = 0; j < nrIntervals; j++) { // SCALE LEVEL
					for (int o = 0; o < in.get(i).get(j).size(); o++) {
						// in.get(i).get(j).get(o).response /= 10;
						Rlist[cnt] = in.get(i).get(j).get(o).response;
						cnt++;
					}
				}
			}
		}
		Arrays.sort(Rlist);
		// System.out.print("SmallesR"+Rlist[100]+" "+Rlist[Rlist.length-1]);
		// return Rlist[1000];
		return Rbarrier;
	}

	public void saveKeyPointImage(ArrayList<ArrayList<ArrayList<KeyPoint>>> in, boolean All, int[] octaves,
			int[] scales) {

		MatOfKeyPoint keypoints_s = new MatOfKeyPoint();
		ArrayList<KeyPoint> temp = new ArrayList<>();
		if (All) {
			for (int i = 0; i < nrOctaves; i++) {// OCTAVE LEVEL
				for (int j = 0; j < nrIntervals; j++) { // SCALE LEVEL
					for (KeyPoint l : in.get(i).get(j)) {
						temp.add(l);
					}
				}
			}
		} else {
			for (int i : octaves) {
				for (int j : scales) {
					for (KeyPoint l : in.get(i).get(j)) {
						temp.add(l);
					}
				}
			}
		}
		keypoints_s.fromList(temp);

		Mat outputImage = new Mat();
		BufferedImage input2 = deepCopy(originalImg);
		// input2 = input2.getSubimage(0, 0, input2.getWidth(),
		// input2.getHeight());

		BufferedImage convertedImg = new BufferedImage(input2.getWidth(), input2.getHeight(),
				BufferedImage.TYPE_3BYTE_BGR);
		convertedImg.getGraphics().drawImage(input2, 0, 0, null);

		byte[] data1 = ((DataBufferByte) convertedImg.getRaster().getDataBuffer()).getData();
		Mat image1 = new Mat(input2.getHeight(), input2.getWidth(), CvType.CV_8UC3);
		image1.put(0, 0, data1);
		Features2d.drawKeypoints(image1, keypoints_s, outputImage);
		String filename = "keypoints.jpg";
		System.out.println(String.format("Writing %s...", filename));
		Highgui.imwrite(filename, outputImage);
		System.out.println(String.format("Wote %s...", filename));

	}

	public BufferedImage makeSmaller(BufferedImage img, int scale) {

		byte[] data = ((DataBufferByte) img.getRaster().getDataBuffer()).getData();
		Mat image = new Mat((int) (img.getHeight()), (int) (img.getWidth()), CvType.CV_8UC3);
		Imgproc.cvtColor(image, image, Imgproc.COLOR_RGB2GRAY);
		image.put(0, 0, data);// variable image gets the data ;
		Mat dst = new Mat(img.getHeight(), img.getWidth(), CvType.CV_8UC3);
		Imgproc.cvtColor(dst, dst, Imgproc.COLOR_RGB2GRAY);

		Size sz = new Size((int) (img.getWidth() / (int) pow(2, scale)), (int) (img.getHeight() / (int) pow(2, scale)));
		Imgproc.resize(image, dst, sz);

		BufferedImage gray = new BufferedImage(dst.width(), dst.height(), BufferedImage.TYPE_BYTE_GRAY);
		byte[] datas = ((DataBufferByte) gray.getRaster().getDataBuffer()).getData();
		dst.get(0, 0, datas);
		return gray;
	}

	public BufferedImage getForPrint() {
		if (forPrint != null) {
			return forPrint;
		}
		return null;
	}

	int n = 0;

	public MatOfKeyPoint getLkeyPoints() {
		MatOfKeyPoint keypoints_s = new MatOfKeyPoint();
		ArrayList<KeyPoint> temp = new ArrayList<>();
		for (int i = 0; i < returnoctave; i++) {// OCTAVE LEVEL
			for (int j = 0; j < returninterval; j++) { // SCALE LEVEL
				for (KeyPoint l : this.openKpL.get(i).get(j)) {
					if ((l.response <= this.Rbarrier || abs(l.response) < 0.5)) {
						n++;
						temp.add(l);
					} else {
						// System.out.println(l.response);
					}
				}
			}
		}
		keypoints_s.fromList(temp);
		System.out.println(n);
		return keypoints_s;
	}

}
