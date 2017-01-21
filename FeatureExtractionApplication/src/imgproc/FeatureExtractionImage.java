package imgproc;

import java.awt.image.BufferedImage;
import java.io.PrintWriter;

import org.opencv.core.Core;
import org.opencv.core.Mat;
import org.opencv.core.MatOfKeyPoint;
import org.opencv.features2d.FeatureDetector;
import org.opencv.features2d.Features2d;
import org.opencv.features2d.KeyPoint;
import org.opencv.highgui.Highgui;

public class FeatureExtractionImage {

	private Mat image;
	private Mat imageKeypoint;

	private MatOfKeyPoint keypoints;

	public FeatureExtractionImage(String imgPath) {
		System.loadLibrary(Core.NATIVE_LIBRARY_NAME);
		image = Highgui.imread(imgPath);
		keypoints = new MatOfKeyPoint();
		imageKeypoint = new MatOfKeyPoint();
		detectKeypointSift();
	}

	private void detectKeypointSift() {
		FeatureDetector detector = FeatureDetector.create(FeatureDetector.SIFT);
		detector.detect(image, keypoints);
		keypoints.toArray();
		KeyPoint[] m = keypoints.toArray();
		try {
			PrintWriter writer = new PrintWriter("the-file-name.txt", "UTF-8");
			writer.print("i.octave");
			writer.print(";");
			writer.print("layer");
			writer.print(";");
			writer.print("octave_extracted");
			writer.print(";");
			writer.print("scale");
			writer.print(";");
			writer.print("size");
			for (KeyPoint i : m) {
				int octave;
				int layer;
				float scale;
				octave = i.octave & 255;
				layer = (i.octave >> 8) & 255;
				octave = octave < 128 ? octave : (-128 | octave);
				scale = octave >= 0 ? 1.f / (1 << octave) : (float) (1 << -octave);
				writer.print(i.octave);
				writer.print(";");
				writer.print(layer);
				writer.print(";");
				writer.print(octave);
				writer.print(";");
				writer.print(scale);
				writer.print(";");
				writer.print(i.size);
			}
			writer.close();
		} catch (Exception e) {

		}
		Features2d.drawKeypoints(image, keypoints, imageKeypoint);
	}

	public Mat getMatImage() {
		return image;
	}

	public Mat getMatImageKeypoint() {
		return imageKeypoint;
	}

	public BufferedImage getBufferedImage() {
		return ImageFormatUtil.convertMatToBufferedImage(image);
	}

	public BufferedImage getBufferedImageKeypoints() {
		return ImageFormatUtil.convertMatToBufferedImage(imageKeypoint);
	}

	public MatOfKeyPoint getKeyPoints() {
		return keypoints;
	}

}
