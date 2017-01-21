package imgproc;

import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;

import org.opencv.calib3d.Calib3d;
import org.opencv.core.Mat;
import org.opencv.core.MatOfPoint2f;
import org.opencv.core.Point;
import org.opencv.core.Rect;
import org.opencv.core.Size;
import org.opencv.imgproc.Imgproc;

public class PanoramaImage {

	FeatureExtractionImage img1, img2;
	FeatureDescriptionImage imgMatch;

	Mat imgResult;

	MatOfPoint2f ptMat1;
	MatOfPoint2f ptMat2;

	public PanoramaImage(FeatureExtractionImage img1, FeatureExtractionImage img2, FeatureDescriptionImage imgMatch) {
		this.img1 = img1;
		this.img2 = img2;
		this.imgMatch = imgMatch;

		ptMat1 = new MatOfPoint2f();
		ptMat2 = new MatOfPoint2f();

		imgResult = new Mat();

		getKeypointsFromGoodMatch();
		computePanoramaImage();
	}

	private void getKeypointsFromGoodMatch() {
		List<Point> ptList1 = new ArrayList<Point>();
		List<Point> ptList2 = new ArrayList<Point>();

		for (int i = 0; i < imgMatch.getGoodMatches().toArray().length; i++) {
			ptList1.add(img1.getKeyPoints().toArray()[imgMatch.getGoodMatches().toArray()[i].queryIdx].pt);
			ptList2.add(img2.getKeyPoints().toArray()[imgMatch.getGoodMatches().toArray()[i].trainIdx].pt);
		}

		ptMat1.fromList(ptList1);
		ptMat2.fromList(ptList2);
	}

	private void computePanoramaImage() {
		// Find the Homography Matrix
		Mat H = Calib3d.findHomography(ptMat2, ptMat1, Calib3d.RANSAC, 1);

		Imgproc.warpPerspective(img2.getMatImage(), imgResult, H,
				new Size(img1.getMatImage().cols() + img2.getMatImage().cols(), img2.getMatImage().rows()));
		Mat half = new Mat(imgResult, new Rect(0, 0, img1.getMatImage().cols(), img1.getMatImage().rows()));
		img1.getMatImage().copyTo(half);

	}

	public Mat getMatPanoramaImage() {
		return imgResult;
	}

	public BufferedImage getBufferedPanoramaImage() {
		return ImageFormatUtil.convertMatToBufferedImage(imgResult);
	}

}
