package imgproc;

import java.util.ArrayList;
import java.util.List;

import org.opencv.core.MatOfKeyPoint;
import org.opencv.features2d.KeyPoint;

public class KeyPointFormatUtil {
	static int octave = 1;
	static int class_id = 1;

	public static MatOfKeyPoint convertArrayToMatOfKeyPoint(ArrayList<double[]> keyPointsOut) {
		MatOfKeyPoint keyPointMat = new MatOfKeyPoint();
		List<KeyPoint> keyPointList = new ArrayList<KeyPoint>();

		KeyPoint keyPoint;
		for (double[] keyPtOut : keyPointsOut) {
			keyPoint = new KeyPoint((float) keyPtOut[0], (float) keyPtOut[1], (float) 2.5, (float) keyPtOut[4],
					(float) keyPtOut[2], octave, class_id);
			keyPointList.add(keyPoint);
		}

		keyPointMat.fromList(keyPointList);

		return keyPointMat;
	}
}
