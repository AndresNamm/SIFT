package imgproc;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import org.opencv.core.Core;
import org.opencv.core.Mat;
import org.opencv.core.MatOfKeyPoint;
import org.opencv.core.Scalar;
import org.opencv.features2d.Features2d;
import org.opencv.highgui.Highgui;

public class FeatureExtractionImage2 extends FeatureExtractionImage {

	private Mat image;
	private Mat imageKeypoint;

	private MatOfKeyPoint keypoints;

	public FeatureExtractionImage2(String imgPath) {
		super(imgPath);
		System.loadLibrary(Core.NATIVE_LIBRARY_NAME);
		image = Highgui.imread(imgPath);
		keypoints = new MatOfKeyPoint();
		imageKeypoint = new MatOfKeyPoint();

		// detectKeypointSift();

		BufferedImage bufImg;
		try {
			bufImg = ImageIO.read(new File(imgPath));
			ImageProcesser imgPro = new ImageProcesser(bufImg, imgPath);
			imgPro.firstPart(bufImg, imgPath);

			// keypoints =
			// KeyPointFormatUtil.convertArrayToMatOfKeyPoint(imgPro.getKeyPointsOut());
			keypoints = imgPro.getLkeyPoints();

			Features2d.drawKeypoints(image, keypoints, imageKeypoint, new Scalar(2, 254, 255),
					Features2d.DRAW_RICH_KEYPOINTS);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

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
